#include "PSO.hpp"

VectorXd ObjectiveFunctionTargetPosition(std::vector<MatrixXd> wrkSpace, Vector3d target, VectorXi* NclosestId)
{
    // Nombre de points à trouver
    const int numClosest = NclosestId->size();
    VectorXd NcloseDist(numClosest);

    std::vector<std::pair<double, int>> distances;

    for (int i = 0; i < wrkSpace.size(); ++i) {
        double dist = (target - wrkSpace[i].rightCols(1)).norm();;
        distances.push_back({dist, i});
    }

    // Trier les distances en ordre croissant (distance la plus faible en premier)
    std::sort(distances.begin(), distances.end(),
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                  return a.first < b.first;
              });
    for (int i = 0; i < numClosest; ++i) {
        (*NclosestId)[i] = distances[i].second;
        NcloseDist[i] = distances[i].first;
    }

    // Retourner la plus petite distance obtenue
    return NcloseDist;
}

VectorXd ObjectiveFunctionObstacle(std::vector<MatrixXd> wrkSpace, VectorXi* NclosestId, Kdtree::KdTree* KdTreeBrain, Target entreeCerveau)
{
    //double minOfMaxs = 99999;
    VectorXd distToDura;
    const int numClosest = NclosestId->size();
    VectorXd NcloseDist(numClosest);
    std::vector<std::pair<double, int>> distances;

    for (int i = 0; i < numClosest; ++i) {
        double distMax = DistMaxCtrDura(wrkSpace[(*NclosestId)[i]], KdTreeBrain, entreeCerveau);
        distances.push_back({distMax, (*NclosestId)[i]});
        //if (distMax < minOfMaxs) minOfMaxs = distMax;
    }

    // Trier les distances en ordre croissant (distance la plus faible en premier)
    std::sort(distances.begin(), distances.end(),
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                  return a.first < b.first;
              });

    for (int i = 0; i < numClosest; ++i) {
        (*NclosestId)[i] = distances[i].second;
        distToDura.conservativeResize(i+1);
        distToDura[i] = distances[i].first;
    } 

    return distToDura;
}

VectorXd ObjectiveFunction(PsoParam PsoParam, std::vector<MatrixXd> wrkSpace, Target target, Kdtree::KdTree* KdTreeBrain, Target entreeCerveau)
{
    VectorXd result(PsoParam.nbObjFunction);
    VectorXi NclosestId(PsoParam.closeTargets);

    VectorXd minDist = ObjectiveFunctionTargetPosition(wrkSpace, target.position, &NclosestId);
    //std::cout << "test ObjectiveFunctionTargetPosition : " << NclosestId[0] << std::endl;
    VectorXd minObst = ObjectiveFunctionObstacle(wrkSpace, &NclosestId, KdTreeBrain, entreeCerveau);

    result << minDist[0], minObst[0];
    return result;
}

std::shared_ptr< geometry::TriangleMesh > PlotValidRobot(PsoParam PSOparam, std::vector<MatrixXd> wrkSpace, std::vector<Vector3d> YtotLastId, std::vector<Target> targets, Kdtree::KdTree* KdTreeBrain, Target entreeCerveau, std::shared_ptr<geometry::TriangleMesh> mesh, Kdtree::KdTree* KdTreeSkull){

    std::shared_ptr<geometry::TriangleMesh> robotsPlot = std::make_shared<open3d::geometry::TriangleMesh>();

    for (int idTarget = 0 ; idTarget < targets.size() ; idTarget++)
    {
        VectorXi NclosestId(1000);
        VectorXd DistToTarget = ObjectiveFunctionTargetPosition(wrkSpace, targets[idTarget].position, &NclosestId);

        if (DistToTarget.size() > 0)
        {
            for (int i=0 ; i < NclosestId.size() ; i++)
            {
                if (DistToTarget[i] > PSOparam.goodDist[0])
                {
                    NclosestId.conservativeResize(i+1);
                    break;
                }
            }
            //std::cout << "NclosestId " << NclosestId << std::endl;

            VectorXd DistToDura = ObjectiveFunctionObstacle(wrkSpace, &NclosestId, KdTreeBrain, entreeCerveau);
            //std::cout << "Target " << idTarget << " dist target : " << DistToTarget[0] << " dist dura : " << DistToDura[0] << std::endl;
            int cpt = 0;
            for (int i=0 ; i < NclosestId.size() ; i++)
            {
                //std::cout << "okayy DistToDura[i] = " << DistToDura[i] << " & PSOparam.goodDistToDura = " << PSOparam.goodDistToDura << std::endl;
                if (DistToTarget.size() > 0 && DistToDura[i] < PSOparam.goodDist[1])
                {
                    //std::cout << "robot isok DistToTarget[i] = " << DistToTarget[i] << " et DistToDura[i] = " << DistToDura[i] << std::endl;
                    //std::cout << "YtotLastId : " << YtotLastId[NclosestId[i]] << std::endl;
                    if (!isRobotContact(wrkSpace[NclosestId[i]], mesh, KdTreeSkull))
                    {
                        std::cout << "Selected for target " << idTarget << " : no " << NclosestId[i] << std::endl;
                        *robotsPlot += *(plotTubes(YtotLastId[NclosestId[i]], wrkSpace[NclosestId[i]]));
                        cpt++;
                    }
                    /**robotsPlot += *(plotTubes(YtotLastId[NclosestId[i]], wrkSpace[NclosestId[i]]));
                    cpt++;*/
                }
            }
            //std::cout << "cible numero " << idTarget << " nb de robots trouvés : " << cpt+1 << std::endl;

            //*robotsPlot += *(plotTubes(VectorXd YtotLastId, MatrixXd posCtr));
        }
        else
            std::cout << "Pas de configuration proche de la cible " << idTarget << std::endl;
        
    }

    return robotsPlot;
}

vector<double> setObjectiveValue(PsoParam PsoParam, vector<double> ctrParam, std::vector<Target> targets, Kdtree::KdTree* KdTreeBrain, Target entreeCerveau,  std::shared_ptr<geometry::TriangleMesh> mesh, Kdtree::KdTree* KdTreeSkull)
{
    std::vector<parameters> vParameters;
    loadParameters("../parameters/parametersMaison.csv",vParameters);
    parameters &pNominal =  vParameters[0];

    std::vector<double> v(PsoParam.nbObjective_values*targets.size(), 99999);

    pNominal.Ux = {ctrParam[0], ctrParam[1], ctrParam[2]};
    pNominal.l_k = {ctrParam[3], ctrParam[4], ctrParam[5]};
    double deltaTubes = 0.005;

    CtrModel ctr(pNominal);

    //WorkspaceData wrkSpace = buildWorkspace(ctr, PsoParam.discretLength, PsoParam.discretAngle, pNominal.offset, deltaTubes, 0.1, false, true, 2*PsoParam.discretAngle);//CtrModel ctr, int discretAngle, int discretLength, VectorXd offset, double deltaTubes, double distUtile, int freqRevolution = 1, bool smart = false
    WorkspaceData wrkSpace = buildWorkspace(ctr, PsoParam.discretAngle, PsoParam.discretLength, pNominal.offset, false, true, 2*PsoParam.discretAngle);
    wrkSpace = cleanWorkspace(wrkSpace, mesh, KdTreeSkull);
    //PlotValidRobot(PsoParam, wrkSpace.data, wrkSpace.YtotLastIdList, targets, KdTreeBrain, entreeCerveau);

    for (int i = 0 ; i < targets.size() ; i++)
    {
        
        VectorXd result = ObjectiveFunction(PsoParam, wrkSpace.data, targets[i], KdTreeBrain, entreeCerveau);
        v[2*i] = result[0];
        v[2*i+1] = result[1];
    }

    return v;
}

// Implémentation de la fonction d'initialisation des particules
void initializeParticles(vector<Particle>& particles, PsoParam param, 
                         const vector<double>& vitesse_max, std::vector<Target> targets, Kdtree::KdTree* KdTreeBrain,
                         Target entreeCerveau, std::shared_ptr<geometry::TriangleMesh> mesh, Kdtree::KdTree* KdTreeSkull) {

    int num_dimensions = param.nbDim;
    std::vector<double> limite_basse = param.lowLimit;
    std::vector<double> limite_haute = param.highLimit;

    for (int i = 0; i < param.nbParticules; ++i) {
        Particle p;
        p.position.resize(num_dimensions);
        p.velocity.resize(num_dimensions);
        p.best_position.resize(num_dimensions);
        p.objective_values.resize(2*targets.size()); // Supposons 2 objectifs
        p.best_objective_values.resize(2*targets.size());

        // Initialisation des positions et vitesses aléatoires
        for (int j = 0; j < num_dimensions; ++j) {
            p.position[j] = limite_basse[j] + (rand() % 1000) / 1000.0 * (limite_haute[j] - limite_basse[j]); // Position entre limite_basse et limite_haute
            p.velocity[j] = (rand() % 200 - 100) / 100.0 * vitesse_max[j]; // Vitesse initiale aléatoire
        }

        std::cout << "ini particule n" << i << " Ux1 = " << p.position[0] << " Ux2 = " << p.position[1] << " Ux3 = " << p.position[2] << " lk1 = " << p.position[3] << " lk2 = " << p.position[4] << " lk3 = " << p.position[5] << std::endl;

        // Calcul des valeurs objectives initiales
        p.objective_values = setObjectiveValue(param, p.position, targets, KdTreeBrain, entreeCerveau, mesh, KdTreeSkull);
        p.best_position = p.position;
        
        particles.push_back(p);
    }
}

// Fonction pour évaluer si une solution est dominée par une autre (Dominance de Pareto)
bool dominates(const std::vector<double>& a, const std::vector<double>& b) {
    bool better = false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] < b[i]) return false;
        if (a[i] > b[i]) better = true;
    }
    return better;
}

// Recherche du meilleur ensemble global basé sur la dominance de Pareto
std::vector<Particle> getParetoFront(std::vector<Particle> oldPareto, std::vector<Particle> particles) {
    std::vector<Particle> candidates;
    candidates.reserve(particles.size() + oldPareto.size()); // optimisation mémoire
    candidates.insert(candidates.end(), particles.begin(), particles.end());
    candidates.insert(candidates.end(), oldPareto.begin(), oldPareto.end());

    //particles.insert(particles.end(), oldPareto.begin(), oldPareto.end());
    vector<Particle> pareto_front;

    for (size_t i = 0; i < candidates.size(); ++i) {
        bool dominated = false;
        for (size_t j = 0; j < candidates.size(); ++j) {
            if (i != j && dominates(candidates[j].objective_values, candidates[i].objective_values)) {
                dominated = true;
                break;
            }
        }
        if (!dominated) {
            pareto_front.push_back(candidates[i]);
        }
    }

    return pareto_front;
}

// Mise à jour des positions et vitesses des particules
void updateParticle(PsoParam param, Particle& p, const vector<double>& global_best_position, 
                    const vector<double>& vitesse_max, std::vector<Target> targets, 
                    Kdtree::KdTree* KdTreeBrain, Target entreeCerveau, 
                    std::shared_ptr<geometry::TriangleMesh> mesh, Kdtree::KdTree* KdTreeSkull) {
    for (size_t i = 0; i < p.position.size(); ++i) {
        // Mise à jour de la vitesse
        p.velocity[i] = 0.5 * p.velocity[i] 
                        + 1.5 * ((rand() % 100) / 100.0) * (p.best_position[i] - p.position[i]) 
                        + 1.5 * ((rand() % 100) / 100.0) * (global_best_position[i] - p.position[i]);
        
        // Limitation de la vitesse
        if (p.velocity[i] > vitesse_max[i]) p.velocity[i] = vitesse_max[i];
        if (p.velocity[i] < -vitesse_max[i]) p.velocity[i] = -vitesse_max[i];

        // Mise à jour de la position
        p.position[i] += p.velocity[i];

        // Limitation de la position dans les limites spécifiées
        if (p.position[i] > param.highLimit[i]) p.position[i] = param.highLimit[i];
        if (p.position[i] < param.lowLimit[i]) p.position[i] = param.lowLimit[i];
        
    }
    std::cout << " Ux1 = " << p.position[0] << " Ux2 = " << p.position[1] << " Ux3 = " << p.position[2] << " lk1 = " << p.position[3] << " lk2 = " << p.position[4] << " lk3 = " << p.position[5] << std::endl;

    // Mise à jour des valeurs objectives après le déplacement
    p.objective_values = setObjectiveValue(param, p.position, targets, KdTreeBrain, entreeCerveau, mesh, KdTreeSkull);

    // Mise à jour de la meilleure position personnelle si nécessaire
    if (dominates(p.objective_values, p.best_objective_values)) {
        p.best_position = p.position;
        p.best_objective_values = p.objective_values;
    }
}

bool isGoodParticle(PsoParam param, int nbTargets, Particle p)
{
    for (int nbFct = 0 ; nbFct < param.nbObjFunction ; nbFct++)
    {
        for (int nbT = 0 ; nbT < nbTargets ; nbT++)
        {
            if (p.objective_values[nbT*param.nbObjFunction + nbFct] > param.goodDist[nbFct])
                return false;
        }
    }
    return true;
}

void updateParticlesWithThreads(std::vector<Particle>& particles,
                                 std::vector<Particle>& result,
                                 const PsoParam& param,
                                 const std::vector<double>& global_best_position,
                                 const std::vector<double>& vitesse_max,
                                 const std::vector<Target>& targets,
                                 Kdtree::KdTree* KdTreeBrain,
                                 Target entreeCerveau,
                                 std::shared_ptr<geometry::TriangleMesh> mesh, 
                                 Kdtree::KdTree* KdTreeSkull) {
    std::mutex result_mutex;
    const size_t num_threads = std::thread::hardware_concurrency();
    const size_t num_particles = particles.size();
    const size_t block_size = (num_particles + num_threads - 1) / num_threads;

    std::vector<std::thread> threads;

    for (size_t t = 0; t < num_threads; ++t) {
        size_t start = t * block_size;
        size_t end = std::min(start + block_size, num_particles);

        threads.emplace_back([&, start, end]() {
            for (size_t i = start; i < end; ++i) {
                updateParticle(param, particles[i], global_best_position, vitesse_max, targets, KdTreeBrain, entreeCerveau, mesh, KdTreeSkull);

                if (isGoodParticle(param, targets.size(), particles[i])) {
                    std::lock_guard<std::mutex> lock(result_mutex);
                    result.push_back(particles[i]);
                }
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }
}

// Fonction principale d'optimisation PSO
vector<Particle> PSO(PsoParam param, std::vector<Target> targets, Kdtree::KdTree* KdTreeBrain, Target entreeCerveau, std::shared_ptr<geometry::TriangleMesh> mesh, Kdtree::KdTree* KdTreeSkull) {
    srand(time(0));

    std::vector<Particle> particles;
    int num_dimensions = param.nbDim;

    // Calcul des vitesses maximales
    vector<double> vitesse_max(num_dimensions);
    for (int i = 0; i < num_dimensions; ++i) {
        vitesse_max[i] = param.highLimit[i] - param.lowLimit[i];
    }

    initializeParticles(particles, param, vitesse_max, targets, KdTreeBrain, entreeCerveau, mesh, KdTreeSkull);

    std::vector<Particle> pareto_front = getParetoFront({}, particles);
    vector<Particle> result = {};


    for (int iter = 0; iter < param.nbIteration; ++iter) {
        
        // Mise à jour des meilleures positions globales
        vector<double> global_best_position(num_dimensions, 0.0);
        if (!pareto_front.empty()) {
            // On prend la meilleure position du front de Pareto (par exemple, la plus proche de l'origine)
            global_best_position = pareto_front[int(pareto_front.size()/2)].position;
        }

        updateParticlesWithThreads(particles, result, param, global_best_position, vitesse_max, targets, KdTreeBrain, entreeCerveau, mesh, KdTreeSkull);

        //-------------sans multithreading
        /*for (size_t i = 0; i < particles.size(); ++i) {
            // Mise à jour des positions et vitesses des particules
            std::cout << "particule n" << i << std::endl;
            updateParticle(param, particles[i], global_best_position, vitesse_max, targets, KdTreeBrain, entreeCerveau);

            if (isGoodParticle(param, targets.size(), particles[i]))
            {
                Particle newPart = particles[i];
                result.push_back(newPart);
            }

        }*/

        //-------------avec multithreading std::execution
        /*std::vector<std::future<std::optional<Particle>>> futures;

        for (size_t i = 0; i < particles.size(); ++i) {
            futures.push_back(std::async(std::launch::async, [&, i]() -> std::optional<Particle> {
                // Mise à jour de la particule
                updateParticle(param, particles[i], global_best_position, vitesse_max, targets, KdTreeBrain, entreeCerveau);

                // Vérification de la validité
                if (isGoodParticle(param, targets.size(), particles[i])) {
                    return particles[i]; // Retourne la particule valide
                } else {
                    return std::nullopt;
                }
            }));
        }

        // Collecte des résultats valides (séquentielle)
        for (auto& fut : futures) {
            auto optPart = fut.get();
            if (optPart.has_value()) {
                result.push_back(optPart.value());
            }
        }*/
        ////-------------fin multithread

        // Recherche du meilleur ensemble Pareto
        vector<Particle> newPareto = getParetoFront(pareto_front, particles);
        pareto_front = newPareto;

        // Affichage des résultats à chaque itération
        std::cout << "Itération " << iter + 1 << ": " << endl;
        int cpt = 1;
        for (const auto& p : pareto_front) {
            std::cout << "pareto front n" << cpt << std::endl;
            std::cout << " Ux1 = " << p.position[0] << " Ux2 = " << p.position[1] << " Ux3 = " << p.position[2] << " lk1 = " << p.position[3] << " lk2 = " << p.position[4] << " lk3 = " << p.position[5] << std::endl;

            std::cout << "obj 1: " << p.objective_values[0] << " | obj 2: " << p.objective_values[1] << " | obj 3: " << p.objective_values[2] << " | obj2 4: " << p.objective_values[3] << std::endl;
            cpt++;
        }
        cout << "-----------------" << endl;
    }
    return result;
}

//-------------------------------------------------------------------------
//-----------------------------  fin PSO  ---------------------------------
//-------------------------------------------------------------------------