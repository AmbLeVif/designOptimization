#include "Workspace.h"

std::vector<Eigen::MatrixXd> genererRevolution(
    const std::vector<Eigen::MatrixXd>& pointsInitiaux,  // Points initiaux
    double discretAngle                                  // Discrétisation de la révolution (en radians)
) {
    std::vector<Eigen::MatrixXd> pointsRevolus;
    // Vérifier que le pas est valide
    double pasAngle = 2*pi/discretAngle;

    // Matrice de rotation autour de l'axe X
    auto rotationX = [](double angle) {
        Eigen::Matrix3d R;
        R << cos(angle), -sin(angle), 0,
             sin(angle),  cos(angle), 0,
             0,           0,          1;
        return R;
    };

    // Itérer sur chaque matrice initiale
    for (const auto& mat : pointsInitiaux) {
        // Vérifier la validité des dimensions (3*n)
        if (mat.rows() != 3) {
            throw std::invalid_argument("Chaque matrice doit avoir 3 lignes.");
        }

        // Nombre de points dans cette matrice
        int n = mat.cols();

        // Ajouter les points initiaux (angle = 0)
        pointsRevolus.push_back(mat);

        // Appliquer les rotations pour chaque angle
        for (double angle = pasAngle; angle < 2*pi; angle += pasAngle) {
            Eigen::Matrix3Xd matRot = rotationX(angle) * mat; // Rotation
            pointsRevolus.push_back(matRot);                 // Ajouter les points tournés
        }
    }
    return pointsRevolus;
}

WorkspaceData buildWorkspace(CtrModel ctr, int discretAngle, int discretLength, VectorXd offset, bool isOriented, bool smart = false, int freqRevolution = -1)
{
    if(freqRevolution == -1)freqRevolution = discretAngle;
    Eigen::Vector<double, NB_TUBES> lTube = ctr.GetTubeParameters().l;
    double l1, l2, l3;
    l1 = lTube[0];
    l2 = lTube[1];
    l3 = lTube[2];

    double deltaP1 = offset[0] - offset[1];
    double deltaP2 = offset[1] - offset[2];
    double deltaP3 = offset[2];

    double deltaD12 = 0.001;
    double deltaD23 = 0.001;
    double deltaD30 = 0.001;

    double eps = 0.00001;

    WorkspaceData workspace;
    std::vector<MatrixXd> result;
    std::vector<Vector3d> YtotLastIdList;
    std::vector<Vector3d> endPos;
    std::vector<MatrixXd> orientations;
    std::vector<VectorXd> listQ;
    VectorXd YtotLastId(NB_TUBES);
    Vector<double,NB_Q> q;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    int cpt = 0;
    int cptError = 0;

    CtrModel ctrSave1 = ctr;
    CtrModel ctrSave2 = ctr;
    CtrModel ctrSave3 = ctr;
    CtrModel ctrSave4 = ctr;
    CtrModel ctrSave5 = ctr;
    int nThread = omp_get_max_threads();
    double alpha = 0;
    for (int ia = 0 ; ia < (smart ? 1 : discretAngle); ia++)
    {
        alpha += 2*pi/discretAngle;
        ctrSave1 = ctr;
        double beta = 0;
        //#pragma omp parallel for num_threads(nThread)

        for(int ib = 0 ; ib < discretAngle ; ib++)
        {
            beta += 2*pi/discretAngle;
            ctrSave2 = ctr;
            double gamma = 0;
            for(int ig = 0 ; ig < discretAngle ; ig++)
            {
                gamma += 2*pi/discretAngle;
                //std::cout << "----------------------------OOOOOH, alpha = " << alpha << ", beta = " << beta << " et gqmma = " << gamma << std::endl;
                ctrSave3 = ctr;

                double start3 = -l3 + deltaD30;
                double end3 = deltaP3;
                int cpt_q3 = 0;
                for (float q3 = start3; q3 <= end3+eps; q3 += (end3-start3)/discretLength)
                {
                    //std::cout << "----------------------------iq3 = " << cpt_q3 << " et q3 = " << q3 << " et " << q3 + (end3-start3)/discretLength << " < " << end3 << std::endl;
                    //std::cout << "----------------------------iq3 = " << cpt_q3 << std::endl;
                    cpt_q3++;
                    ctrSave4 = ctr;

                    double start2 = q3 + l3 + deltaD23 -l2 + eps; 
                    double end2 = q3 + deltaP2;
                    int cpt_q2 = 0;
                    for (float q2 = start2; q2 <= end2+eps; q2 += (end2-start2)/discretLength)
                    {
                        //std::cout << "----------------------------iq2 = " << cpt_q2 << std::endl;
                        //std::cout << "----------------------------iq2 = " << cpt_q2 << " et q2 = " << q2 << " et " << q2 + (end2-start2)/discretLength << " < " << end2 << std::endl;
                        cpt_q2++;
                        ctrSave5 = ctr;
                        double start1 = q2 + l2 + deltaD12 -l1 + eps; 
                        double end1 = q2 + deltaP1;
                        int cpt_q1 = 0;
                        for (float q1 = start1; q1 <= end1+eps; q1 += (end1-start1)/discretLength)
                        {
                            //std::cout << "----------------------------iq1 = " << cpt_q1 << std::endl;
                            cpt_q1++;
                            q << q1, q2, q3, gamma, beta, alpha; // coord. articulaires (angle en rad) pi/4
                            // Compute model
                            //std::cout << "ah" << std::endl;
                            if (ctr.Compute(q, opt_LOAD) >=0)
                            {
                                workspace.isStable.push_back(1);
                                for (int i = 0; i<NB_TUBES; i++){YtotLastId[i]=getYtotIndexFromIend(ctr.segmented.iEnd(i)) + 1;}
                                //std::cout << "wiiii cols = " << ctr.GetYTot().cols() << " rows = " << ctr.GetYTot().rows() << " enfin YtotLastId= " << YtotLastId << std::endl;

                                MatrixXd formeCtr = ctr.GetYTot().topLeftCorner(3, YtotLastId[0]);
                                result.push_back(formeCtr);
                                YtotLastIdList.push_back(YtotLastId);
                                listQ.push_back(q);
                                if (isOriented)
                                {
                                    Eigen::Matrix3d orientationMatrix;
                                    for (int i = 0; i < NB_TUBES; ++i) {
                                        orientationMatrix(0, i) = ctr.GetYTot()(i + 3, getYtotIndexFromIend(ctr.segmented.iEnd(0)));    // Ligne 3 à 5 (index 2 à 4)
                                        orientationMatrix(1, i) = ctr.GetYTot()(i + 6, getYtotIndexFromIend(ctr.segmented.iEnd(0)));    // Ligne 6 à 8 (index 5 à 7)
                                        orientationMatrix(2, i) = ctr.GetYTot()(i + 9, getYtotIndexFromIend(ctr.segmented.iEnd(0)));    // Ligne 9 à 11 (index 8 à 10)
                                        orientations.push_back(orientationMatrix);
                                    }
                                }
                                else
                                    orientations.push_back(Matrix3d::Identity());
                            }
                            else
                            {
                                cptError++;
                                workspace.isStable.push_back(0);
                            }
                            cpt++;
                        }
                        ctr = ctrSave5;

                    }
                    ctr = ctrSave4;
                }
                ctr = ctrSave3;
            }
            ctr = ctrSave2;
        }
        ctr = ctrSave1;
    }
    result.erase(result.begin());
    YtotLastIdList.erase(YtotLastIdList.begin());
    listQ.erase(listQ.begin());

    std::vector<Vector3d> YtotLastIdListRevolution;
    std::vector<int> isStableRevolution;
    std::vector<VectorXd> listQRevolution;
    if (freqRevolution > 1 && smart)
    {
        result = genererRevolution(result, freqRevolution);
        orientations = genererRevolution(orientations, freqRevolution);
        
        for (auto value : YtotLastIdList) {for (int i = 0; i < freqRevolution; ++i) {YtotLastIdListRevolution.push_back(value);} }
        for (const MatrixXd& mat : result){endPos.push_back(mat.col(mat.cols() - 1));}
        for (auto value : workspace.isStable) {for (int i = 0; i < freqRevolution; ++i) {isStableRevolution.push_back(value);} }
        for (auto value : listQ) {
            for (int i = 0; i < freqRevolution; ++i) {
                value[3] = 2*pi/freqRevolution*(i+1);
                value[4] = value[4]+pi;
                value[5] = value[5]+pi;
                
                listQRevolution.push_back(value);
            } 
        }
        
        YtotLastIdList = YtotLastIdListRevolution;
        workspace.isStable = isStableRevolution;
        listQ = listQRevolution;
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()*1e-6 << "[s] pour cpt = " << cpt << std::endl;
    std::cout << "size data " << result.size()  << "size stable " << isStableRevolution.size() << "lastId stable " << YtotLastIdList.size() << std::endl;

    workspace.q = listQ;
    workspace.data = result;
    workspace.orientations = orientations;
    workspace.YtotLastIdList = YtotLastIdList;
    workspace.endPos = endPos;
    return workspace;
}

int main_wrkspace(int argc, char *argv[]) {

    //Test robot ------------------------------
    std::vector<parameters> vParameters;
    if(loadParameters("../parameters/parameters.csv",vParameters) != 0){
      return -1;
    }

    parameters &pNominal =  vParameters[0];

    double l1 = pNominal.l[0];
    double l2 = pNominal.l[1];
    double l3 = pNominal.l[2];
    std::cout << "l1 : " << pNominal.l[0] << " l2 : " << pNominal.l[1] << " l3 : " << pNominal.l[2] << std::endl;

    CtrModel ctr(pNominal);


    Vector<double,NB_Q> q;
    
    q << -l1+0.003,-l2+0.002,-l3+0.001, 0, pi, 0;
    // Compute model
    ctr.Compute(q, opt_LOAD);


    int discretAngle = 2;
    int discretLength = 10;

    //CtrModel ctr, int discretAngle, int discretLength, VectorXd offset, bool isOriented, bool smart = false, int freqRevolution = 1
    WorkspaceData wrkSpace = buildWorkspace(ctr, discretAngle, discretLength, pNominal.offset,false, true);

    std::cout << "taille wrksp : " << wrkSpace.data.size() << std::endl;

    return 0;
}