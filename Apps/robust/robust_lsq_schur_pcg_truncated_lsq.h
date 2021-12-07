#ifndef ROBUST_LSQ_SCHUR_PCG_TRUNC_LSQ
#define ROBUST_LSQ_SCHUR_PCG_TRUNC_LSQ

#include "robust_lsq_common.h"
#include "Math/v3d_linearbase.h"
#include "utilities_common.h"
#include "Math/v3d_linear_lu.h"
#include "Math/v3d_linear_ldlt.h"

#include "robust_lsq_schur_pcg_common.h"

#include <iomanip>
#include <fstream>
#include <random>

namespace Robust_LSQ{
    struct Robust_LSQ_Optimizer_Schur_PCG_TruncLSQ : public Robust_LSQ_Optimizer_Schur_PCG_Base{
        typedef Robust_LSQ_Optimizer_Schur_PCG_Base Base;

        //Constructor
        Robust_LSQ_Optimizer_Schur_PCG_TruncLSQ(NLSQ_ParamDesc const &paramDesc, 
                std::vector<NLSQ_CostFunction *> const &costFunctions,
                std::vector<Robust_NLSQ_CostFunction_Base *> const &robustCostFunctions):
            Base(paramDesc, costFunctions, robustCostFunctions)
        {

        }


        //Compute the robust cost given a vector of errors
        double compute_robust_cost (Vector<double> const &errors) {
            double cost = 0;
            for (int i=0; i < _nMeasurements; i++)
                cost += _robustCostFunctions[0]->eval_target_fun(1.0, errors[i]);
            /* cost +=  (errors[i] > 1.0 ? 1 : errors[i]); */
            return cost;
        }

        double compute_weighted_cost (Vector<double> const &errors, Vector<double> const &weights) {
            double cost = 0;
            for (int i=0; i < _nMeasurements; i++){
                double const x = errors[i];
                /* cost += weights[i] * errors[i]; */
                /* cost+= (abs(x) > 1.0 ? 1 : x); */
                cost += clamp(x);
            }
            return cost;
        }



        double clamp(double const x) const {
            if (x > 1.0)
                return (1.0 - _slope) + _slope * x;
            else if (x >=0 || x<=1 )
                return x;
            else 
                return 0 +  _slope * x ;
            /* return (abs(x) > 1.0 ? 1 + _slope * x : x); */
            /* return (abs(x) > 1.0 ? 1 + _slope * x : x); */
        }

        double hardtanh (double const x) const {
            if (x > 1.0) return (1.0 - _slope) + _slope * x;
            else if (x < -1.0) return -1.0;
            else return x;
        }

        void minimize(){
            int const totalParamDimension = _totalParamDimension;
            int const nObjs = _robustCostFunctions.size();
            double &damping_value = this->_damping ;

            vector<Vector<double>> cached_errors(nObjs);
            vector<Vector<double>> cached_weights(nObjs);
            for (int obj = 0; obj < nObjs; ++obj) {
                cached_errors[obj].newsize(_nMeasurements);
                cached_weights[obj].newsize(_nMeasurements);
                fillVector(1.0, _residuals[obj]->_weights);
            }

            for (int obj = 0; obj < nObjs; ++obj)
                _robustCostFunctions[obj]->cache_residuals(
                        *_residuals[obj], cached_errors[obj]);

            _timer.reset(); _timer.start();
            _bestCost = 1e20;


            /* std::vector<double> scale_list{50, 25, 15, 5, 1}; */
            /* std::vector<double> scale_list{5, 2, 1}; */
            std::vector<double> slope_list{0.00};
            for (int outer_iter = 0; outer_iter < slope_list.size(); outer_iter++){
                /* _scale = scale_list[outer_iter]; */
                _slope=slope_list[outer_iter];

                for (int iter = 0; iter < maxIterations; iter++){

                    if (iter % 1==0) {
                        /* _tau *= 0.8; */
                        _mu *= 1.05;
                        /* _slope *= 0.8; */
                    }

                    //Cache residuals 
                    _robustCostFunctions[0]->cache_residuals(
                            *_residuals[0], cached_errors[0]);

                    //Compute robust cost
                    double const initCost = compute_robust_cost(cached_errors[0]);
                    double const initTruncatedCost = computeTruncatedLSQ(cached_errors[0]);
                    this->updateStat(initCost, initTruncatedCost);

                    //Now, we use our special way to compute weights
                    for (int k=0; k < _nMeasurements; k++){
                        double ek = cached_errors[0][k];
                        double const zstar = clamp(ek);
                        cached_weights[0][k] = 0.5;
                        /* double tauk = 0.001; */
                        if (ek > 1.0){
                            double tauk = 1/_mu*(2*(ek -1) + 4 * (ek-1)*(ek-1)*_mu * _mu + 2 * _mu * (ek-1));
                            if (tauk <= ek -1)
                                tauk = 1.0;

                            double const zbar = clamp((ek + _mu * tauk * zstar - 0.5 * tauk)/(1 + _mu * tauk));
                            cached_weights[0][k] = sqrt((zstar - zbar)/tauk);
                        }

                    }

                    std::cout << "iter : " << iter << 
                        " TLS Cost = " << initTruncatedCost << 
                        " damping = " << damping_value << endl;

                    this->fillJacobians();
                    this->evalJt_e(cached_weights);
                    this->fillHessian(cached_weights);
                    this->addDamping();
                    this->prepare_schur();
                    this->fillJt_e_schur();
                    this->fill_M();
                    this->solveJtJ();

                    //Update parameters
                    this->saveAllParameters();
                    if (isnan(_deltaA[0][0])){
                        damping_value *= 10;
                        continue;
                    }
                    this->updateParameters(0, _deltaA);
                    this->updateParameters(1, _deltaB);
                    this->finishUpdateParameters();

                    _robustCostFunctions[0]->cache_residuals(*_residuals[0], cached_errors[0]);
                    double const newCost = compute_robust_cost(cached_errors[0]);

                    double const newTruncatedCost = computeTruncatedLSQ(cached_errors[0]);
                    bool accepted = false, converged = false;

                    if (newTruncatedCost <= initTruncatedCost){
                        /* if (newCost <= initCost){ */
                        accepted = true;
                        damping_value = max(1e-8, damping_value * 0.1);
                        if (abs(newCost - initCost) < 0.01 * initCost){
                            converged = true;
                        }
                    } else {
                        if (abs(initTruncatedCost - newTruncatedCost) < 1e-9)
                            break;
                        damping_value *= 10;
                        this->restoreAllParameters();
                        iter--;
                    }
                    } //for iter
                }

            }

            private:
            double _tau = 5.0;
            double _slope = 0.00;
            double _scale = 0.5; // Potentially used for GNC-type algorithms 
            double _mu = 1.0;

        };
    }// end namespace Robust_LSQ




#endif 

