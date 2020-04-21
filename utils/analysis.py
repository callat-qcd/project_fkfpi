#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np


class BayesModelAvg:
    ''' Under Bayes Theorem, models which are used to fit the same dataset
        can easily be averaged together with the exp(logGBF) factor providing
        a weight, such that the relative weigth of any given model with respect
        to the most likely model (of those selected) is given by

        w_i = exp( logGBF_i - logGBF_max )

        The final uncertainty under such model averaging is given simply by a
        combination of the weighted variance from each fit, plus the variance
        over the models:

        mu  = sum_i w_i mu_i
        Var = sum_i w_i * var_i   + [ sum_j mu_j**2 * w_j - mu**2 ]

        This class performs this weighted average, and creates histograms of
        the final answer including splitting the historgram in various ways to
        understand the dominant contributions to the model average uncertainty
    '''
    def __init__(self, results):
        ''' results is a dictionary where the key is the model name and the
            value is an lsqfit object
        '''
        self.results = results
        # create a list of the keys so we gaurantee the same order
        self.r_list  = [a_res for a_res in self.results]
        self.weights = self.get_weights()

    def get_weights(self):
        weights = []
        for a_res in self.r_list:
            weights.append(np.exp(self.results[a_res].logGBF))
        weights = np.array(weights)
        weights = weights / np.max(weights)
        weights = weights / weights.sum()
        return weights

    def print_weighted_models(self):
        i_weights = np.argsort(self.weights)[::-1]
        print(r"model & $\chi^2/{\rm dof}$ & $Q$& logGBF& weight& $F_K/F_\pi$\\")
        print(r'\hline')
        for a_i, a_model in enumerate(np.array(self.r_list)[i_weights]):
            chi2   = self.results[a_model].chi2
            dof    = self.results[a_model].dof
            Q      = self.results[a_model].Q
            logGBF = self.results[a_model].logGBF
            w_i    = self.weights[i_weights][a_i]
            phys   = self.results[a_model].phys
            print(r'%33s &  %.3f&  %.3f&  %.3f&  %.3f&  %s\\' %(a_model, chi2/dof, Q, logGBF, w_i, phys))
