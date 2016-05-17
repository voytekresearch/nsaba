"""
visualizer.py:
Visualization tools for nsaba module.
Author: Torben Noto & Simon Haxby
"""
from nsaba import Nsaba
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from nsabatools import not_operational, preprint


class NsabaVisualizer(object):
    """
    Wraps around Nsaba object, providing visualization methods.
    """
    def __init__(self, nsaba_obj):
        if type(nsaba_obj) == Nsaba:
            self.no = nsaba_obj
        else:
            raise ValueError("NsabaVisualizer() parameter not a Nsaba instance")
    
    def visualize_ge(self, gene, alpha=0.4, figsize=(16, 10), z_score=True):
        """
        Generates 3-D heat map of gene expression for a specified gene.

        Parameters
        ----------
        gene : str
            Entrez ID of gene for heat map generation.

        alpha : float
            Sets color density of coordinates used in heat map.
        """
        if gene in self.no.ge:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111, projection='3d')
            if z_score is True:
                weights = (np.nan_to_num(self.no.ge[gene]["mean"]['GE'])-np.nanmean(self.no.ge[gene]["mean"]['GE']))/np.nanstd(self.no.ge[gene]["mean"]['GE'])
            else:
                weights = np.nan_to_num(self.no.ge[gene]["mean"]['GE'])
            colors = cm.viridis(weights)
            color_map = cm.ScalarMappable(cmap=cm.viridis)
            color_map.set_array(weights)
            fig.colorbar(color_map)

            x = self.no._aba['mni_coords'].data[:, 0]
            y = self.no._aba['mni_coords'].data[:, 1]
            z = self.no._aba['mni_coords'].data[:, 2]

            ax.scatter(x, y, z, c=colors, alpha=alpha)
            ax.set_title('Gene Expression of gene ID ' + str(gene))
            return fig
        else:
            raise ValueError("Gene %s has not been initialized. "
                             "Use self.no.get_aba_ge([%s])" % str(e))

    def visualize_ns(self, term, alpha=0.5, figsize=(16, 10), z_score=True):

        if term in self.no.term:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111, projection='3d')
            if z_score is True:
                weights = (np.nan_to_num(no.term[term]['act'])-np.nanmean(self.no.term[term]['act']))/np.nanstd(self.no.term[term]['act'])
            else:
                weights = np.nan_to_num(self.no.term[term]['act'])
            colors = cm.viridis(weights)
            color_map = cm.ScalarMappable(cmap=cm.viridis)
            color_map.set_array(weights)
            fig.colorbar(color_map)

            x = self.no._aba['mni_coords'].data[:, 0]
            y = self.no._aba['mni_coords'].data[:, 1]
            z = self.no._aba['mni_coords'].data[:, 2]

            ax.scatter(x, y, z, c=colors, alpha=alpha)
            ax.set_title('Map of ' + term)
            return fig
        else:
            raise ValueError("Term %s has not been initialized. "
                             "Use N.get_ns_act([%s])" % term)

    def visualize_ns_all(self, term, no_ids=10, alpha=0.2, figsize=(16, 10), z_score=True):
        """
        Generates 3-D heat map of term activation for a specified term.

        Parameters
        ----------
        gene : str
            Entrez ID of gene for heat map generation.

        no_ids : int
            Sets specificity of heat map. // More detail needed

        alpha : float
            Sets color density of coordinates used in heat map.
        """
        if term in self.no.term:
            try:
                len(self.no._ns['id_dict'])
            except KeyError:
                self.no.ns_load_id_dict()
            heat = (self.no._ns['features_df'][term]-self.no._ns['features_df'][term].mean())/self.no._ns['features_df'][term].std()
            sorted_heat_vals = sorted(enumerate(heat), key=lambda x: x[1], reverse=True)[0:no_ids]
            weights = zip(*sorted_heat_vals)[1]
            inds = zip(*sorted_heat_vals)[0]
            pmids = [self.no._ns['features_df']['pmid'].ix[ind] for ind in inds]
            all_coords = []
            for pmid in pmids:
                if len(self.no._ns['database_df'].loc[self.no._ns['database_df']['id'] == pmid]) > 0:
                    all_coords.append(self.no._ns['id_dict'][pmid])
            xvals = []
            yvals = []
            zvals = []
            new_weights = []
            wc = 0
            for coord_set in all_coords:
                wc += 1
                for coord in coord_set:
                    xvals.append(coord[0])
                    yvals.append(coord[1])
                    zvals.append(coord[2])
                    new_weights.append(weights[wc])
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111, projection='3d')
            colors = cm.viridis(new_weights/max(new_weights))
            color_map = cm.ScalarMappable(cmap=cm.viridis)
            color_map.set_array(new_weights)
            fig.colorbar(color_map)
            ax.scatter(xvals, yvals, zvals, c=colors, alpha=alpha)
            ax.set_title('Heat map of ' + str(term))
            return fig
        else:
            raise ValueError("Term '%s' has not been initialized. Use get_ns_act('%s')" % term)

    @not_operational
    def visualize_ns_old(self, term, points=200):
        """
        Use randomly selected coordinates instead of most active
        """
        if term in self.no.term:
            term_index = self.no._ns['features_df'].columns.get_loc(term)
            rand_point_inds = np.random.random_integers(0, len(np.squeeze(zip(self.no._ns['mni_coords'].data))), points)
            rand_points = np.squeeze(zip(self.no._ns['mni_coords'].data))[rand_point_inds]
            weights = []
            inds_of_real_points_with_no_fucking_missing_study_ids = []
            for rand_point in range(len(rand_points)):
                if len(self.no.coord_to_ns_act(rand_points[rand_point].astype(list))) > 0:
                    inds_of_real_points_with_no_fucking_missing_study_ids.append(rand_point_inds[rand_point])
                    weights.append(self.no.coord_to_ns_act(rand_points[rand_point].astype(list))[term_index])
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            colors = cm.viridis(weights/max(weights))
            color_map = cm.ScalarMappable(cmap=cm.viridis)
            color_map.set_array(weights)
            fig.colorbar(color_map)
            x = self.no._ns['mni_coords'].data[inds_of_real_points_with_no_fucking_missing_study_ids, 0]
            y = self.no._ns['mni_coords'].data[inds_of_real_points_with_no_fucking_missing_study_ids, 1]
            z = self.no._ns['mni_coords'].data[inds_of_real_points_with_no_fucking_missing_study_ids, 2]
        else:
            raise ValueError('Term '+term + ' has not been initialized. '
                                            'Use get_ns_act(' + term + ')')
        ax.scatter(x, y, z, c=colors, alpha=0.4)
        ax.set_title('Estimation of ' + term)

    @not_operational
    def lstsq_ge_ns(self, gene, term, logy=False, logx=False, regression_line=True):
        self.lstsq_ns_ge(self, term, gene, logy=logy, logx=logx, regression_line=True)

    def lstsq_ns_ge(self, term, gene, logy=False, logx=False, only_term=False, verbose=False, regression_line=True):
        """
        Generates a linear regression plot for gene expression to
        term activation.

        Least squares solution is used.

        Parameters
        ----------
        term : str
            NS term to be linearly correlated/predicted with a given gene's
            gene expression.

        gene : str
            Entrez ID of gene whose gene expression is to be used as
            a linear predictor of term activation.

        logy : bool
            log-space of y-axis.

        logx : bool
            log-space of x-axis.

         Returns
        -------
        correlation, [m,c] : np.array([2 x 2]), [2]
            Returns a tuple consisting of the gene-term correlation matrix
            and regression parameters.

        """
        for g in gene:
            if g in self.no.ge:
                if term in self.no.term:
                    ge_ns_mat = self.no.matrix_builder([term], gene)
                    if only_term:
                        if ge_ns_mat.shape[0] > 900:  # check this num later
                            print 'reinitializing ' + term + ' for hacky plotting method'
                            self.no.est_ns_act(term, radius=10)
                            ge_ns_mat = self.no.matrix_builder([term], gene)
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    ax.grid(False)
                    if logy:
                        ax.set_yscale('log')
                    if logx:
                        ax.set_xscale('log')

                    # remove nans
                    valid1 = np.isfinite(ge_ns_mat[:, 0])
                    valid2 = np.isfinite(ge_ns_mat[:, 1])
                    valid_inds = np.logical_and(valid1, valid2)

                    # correlation
                    r, p = np.corrcoef(ge_ns_mat[:, 0][valid_inds], ge_ns_mat[:, 1][valid_inds])

                    # linear regression
                    X = np.vstack([ge_ns_mat[:, 0][valid_inds], np.ones(len(ge_ns_mat[:, 0][valid_inds]))]).T
                    m, c = np.linalg.lstsq(X, ge_ns_mat[:, 1][valid_inds])[0]

                    # print 'Correlation between ' + term + ' and gene number ' + str(gene)
                    # print correlation
                    # print 'Linear regression between ' + term + ' and gene number ' + str(gene)
                    #       +' Slope =' + str(m) + ' y intercept = '+ str(c)

                    ax.plot(ge_ns_mat[:, 0], ge_ns_mat[:, 1], '.')
                    ax.grid(False)
                    if regression_line is True:
                        ax.plot([min(ge_ns_mat[:, 0]), max(ge_ns_mat[:, 0][valid_inds])],
                                [m*min(ge_ns_mat[:, 0][valid_inds])+c,
                                m*max(ge_ns_mat[:, 0][valid_inds])+c], 'r')
                    ax.set_xlabel(str(gene))
                    ax.set_ylabel(term)
                    if verbose is True:
                        print 'correlation: r=' + str(r[1]) + '\n'
                        print 'linear regression: m=' + str(m) + ' c=' + str(c)
                    return [r, p], [m, c]
                else:
                    raise ValueError("Term '%s' has not been initialized. Use get_ns_act('%s')" % term)
            else:
                raise ValueError("Gene %s has not been initialized. "
                                 "Use self.no.get_aba_ge([%s])" % str(g))

    def lstsq_ns_ns(self, term1, term2, logy=False, logx=False):
        """
        Generates a linear regression plot for term activation to
        term activation.

        Least squares solution is used.

        Parameters
        ----------
        term1 : str
            NS predictor term.

        term2 : str
            NS term to be linearly correlated/predicted with a predictor term's
            term activation.

        logy : bool
            log-space of y-axis.

        logx : bool
            log-space of x-axis.

        Returns
        -------
        correlation, [m,c] : np.array([2 x 2]), [2]
            Returns a tuple consisting of the term-term correlation matrix
            and regression parameters.

        """
        if term1 in self.no.term:
            if term2 in self.no.term:
                # remove nans
                valid1 = np.isfinite(self.no.term[term1]['act'])
                valid2 = np.isfinite(self.no.term[term1]['act'])
                valid_inds = np.logical_and(valid1, valid2)

                # Correlation
                correlation = np.corrcoef(self.no.term[term1]['act'][valid_inds], self.no.term[term2]['act'][valid_inds])
                # linear Regression
                regression_matrix = np.vstack([self.no.term[term1]['act'][valid_inds],
                                               np.ones(len(self.no.term[term1]['act'][valid_inds]))]).T
                m, c = np.linalg.lstsq(regression_matrix, self.no.term[term2]['act'][valid_inds])[0]

                # Plotting
                fig = plt.figure()
                ax = fig.add_subplot(111)
                if logy:
                    ax.set_yscale('log')
                if logx:
                    ax.set_xscale('log')
                ax.plot(self.no.term[term1]['act'][valid_inds], self.no.term[term2]['act'][valid_inds], '.')
                ax.plot([min(self.no.term[term1]['act'][valid_inds]), max(self.no.term[term1]['act'][valid_inds])],
                        [m*min(self.no.term[term2]['act'][valid_inds])+c,
                         m*max(self.no.term[term2]['act'][valid_inds])+c], 'r')
                ax.set_xlabel(term1)
                ax.set_ylabel(term2)
                return correlation, [m, c]
            else:
                raise ValueError("Term '%s' has not been initialized. Use get_ns_act('%s')" % term2)
        else:
            raise ValueError("Term '%s' has not been initialized. Use get_ns_act('%s')" % term1)

    def lstsq_ge_ge(self, gene1, gene2):
        """
        Generates a linear regression plot for gene expression
        to gene expression.

        Least squares solution is used.

        Parameters
        ----------
        gene1 : str
            Entrez ID of ABA predictor gene.

        gene2 : str
            Entrez ID of ABA gene to be linearly correlated/predicted with a predictor gene's
            gene expression.

        logy : bool
            log-space of y-axis.

        logx : bool
            log-space of x-axis.

        Returns
        -------
        correlation, [m,c] : np.array([2 x 2]), [2]
            Returns a tuple consisting of the gene-gene correlation matrix
            and regression parameters.

        """

        genes = [gene1, gene2]

        for gene in genes:
            if gene not in self.no.ge:
                raise ValueError("Gene %s has not been initialized. "
                                 "Use self.no.get_aba_ge([%s])" % str(gene))

        if len(self.no.ge[genes[0]]["mean"]['GE']) != len(self.no.ge[genes[1]]["mean"]['GE']):
            raise ValueError("'GE' size mismatched rerun Nsaba.estimate_aba_ge() again.")

        g0 = self.no.ge[genes[0]]["mean"]['GE']
        g1 = self.no.ge[genes[1]]["mean"]['GE']
        # Correlation
        correlation = np.corrcoef(g0, g1)
        # linear regression
        regression_matrix = np.vstack([g0, np.ones(len(g0))]).T
        m, c = np.linalg.lstsq(regression_matrix, g1)[0]

        # Plotting
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.plot(g0, g1, '.')
        plt.plot([min(g0), max(g0)],
                 [m*min(g1)+c, m*max(g1)+c], 'r')
        ax.set_xlabel(genes[0])
        ax.set_ylabel(genes[1])
        return correlation, [m, c]
