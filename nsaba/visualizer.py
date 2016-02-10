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
    
    def visualize_ge(self, gene, alpha=0.4):
        """
        Generates 3-D heat map of gene expression for a specified gene.

        Parameters
        ----------
        gene : str
            Entrez ID of gene for heat map generation.

        alpha : float
            Sets color density of coordinates used in heat map.
        """
        for e in gene:
            if e in self.no.ge:
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                weights = self.no.ge[e]["mean"]['GE']
                colors = cm.jet(weights/max(weights))
                color_map = cm.ScalarMappable(cmap=cm.jet)
                color_map.set_array(weights)
                fig.colorbar(color_map)

                x = self.no._aba['mni_coords'].data[:, 0]
                y = self.no._aba['mni_coords'].data[:, 1]
                z = self.no._aba['mni_coords'].data[:, 2]

                ax.scatter(x, y, z, c=colors, alpha=alpha)
            else:
                raise ValueError("Gene %s has not been initialized. Use self.no.get_aba_ge([%s])" % (str(e), str(e)))
        ax.set_title('Gene Expression of gene ID ' + str(gene))

    def visualize_ns(self, term, no_ids=10, alpha=0.2):
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
            heat = self.no._ns['features_df'][term]
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
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            colors = cm.jet(new_weights/max(new_weights))
            color_map = cm.ScalarMappable(cmap=cm.jet)
            color_map.set_array(new_weights)
            fig.colorbar(color_map)
            ax.scatter(xvals, yvals, zvals, c=colors, alpha=alpha)
            ax.set_title('Heat map of ' + str(term))
        else:
            raise ValueError("Term '%s' has not been initialized. Use get_ns_act('%s')" % (term,term))

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
            colors = cm.jet(weights/max(weights))
            color_map = cm.ScalarMappable(cmap=cm.jet)
            color_map.set_array(weights)
            fig.colorbar(color_map)
            x = self.no._ns['mni_coords'].data[inds_of_real_points_with_no_fucking_missing_study_ids, 0]
            y = self.no._ns['mni_coords'].data[inds_of_real_points_with_no_fucking_missing_study_ids, 1]
            z = self.no._ns['mni_coords'].data[inds_of_real_points_with_no_fucking_missing_study_ids, 2]
        else:
            raise ValueError('Term '+term + ' has not been initialized. Use get_ns_act(' + term + ')')
        ax.scatter(x, y, z, c=colors, alpha=0.4)
        ax.set_title('Estimation of ' + term)

    def lstsq_ge_ns(self, gene, term, logy=False, logx=False):
        """
        Calls lstsq_ns_ge
        """
        self.lstsq_ns_ge(self, term, gene, logy=logy, logx=logx)


    def lstsq_ns_ge(self, term, gene, logy=False, logx=False, only_term=False, verbose=False):
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
                    if logy:
                        ax.set_yscale('log')
                    if logx:
                        ax.set_xscale('log')

                    # correlation
                    correlation = np.corrcoef(ge_ns_mat[:, 0], ge_ns_mat[:, 1])

                    # linear regression
                    X = np.vstack([ge_ns_mat[:, 0], np.ones(len(ge_ns_mat[:, 0]))]).T
                    m, c = np.linalg.lstsq(X, ge_ns_mat[:, 1])[0]

                    ax.plot(ge_ns_mat[:, 0], ge_ns_mat[:, 1], '.')
                    ax.plot([min(ge_ns_mat[:, 0]), max(ge_ns_mat[:, 0])], [m*min(ge_ns_mat[:, 0])+c,
                                                                           m*max(ge_ns_mat[:, 0])+c], 'r')
                    ax.set_xlabel(str(gene))
                    ax.set_ylabel(term)
                    if verbose:
                        print 'Correlation: r=' + str(correlation[0])+' p=' + str(correlation[1]) + \
                              '\n Linear regression coefficients: ' + 'm=' + str(m) + ' c=' + str(c)

                    return correlation, [m, c]
                else:
                    raise ValueError("Term '%s' has not been initialized. Use get_ns_act('%s')" % (term, term))
            else:
                raise ValueError("Gene %s has not been initialized. Use self.no.get_aba_ge([%s])" % (str(g),str(g)))

    def lstsq_ns_ns(self, term1, term2, logy=False, logx=False, verbose=False):
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
                # Correlation
                correlation = np.corrcoef(self.no.term[term1]['act'], self.no.term[term2]['act'])
                # linear Regression
                regression_matrix = np.vstack([self.no.term[term1]['act'],
                                               np.ones(len(self.no.term[term1]['act']))]).T
                m, c = np.linalg.lstsq(regression_matrix, self.no.term[term2]['act'])[0]

                # Plotting
                fig = plt.figure()
                ax = fig.add_subplot(111)
                if logy:
                    ax.set_yscale('log')
                if logx:
                    ax.set_xscale('log')
                ax.plot(self.no.term[term1]['act'], self.no.term[term2]['act'], '.')

                #linear regression line
                ax.plot([min(self.no.term[term1]['act']), max(self.no.term[term1]['act'])],
                        [m*min(self.no.term[term2]['act'])+c,
                         m*max(self.no.term[term2]['act'])+c], 'r')
                ax.set_xlabel(term1)
                ax.set_ylabel(term2)
                if verbose:
                    print 'Correlation: r=' + str(correlation[0])+' p=' + str(correlation[1]) + \
                                '\n Linear regression coefficients: ' + 'm=' + str(m) + ' c=' + str(c)
                return correlation, [m, c]
            else:
                raise ValueError("Term '%s' has not been initialized. Use get_ns_act('%s')" % (term2, term2))
        else:
            raise ValueError("Term '%s' has not been initialized. Use get_ns_act('%s')" % (term1, term1))

    def lstsq_ge_ge(self, gene1, gene2, verbose=False):
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
                raise ValueError("Gene %s has not been initialized. Use self.no.get_aba_ge([%s])" % (str(gene), str(gene)))

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
        if verbose:
            print 'Correlation: r=' + str(correlation[0])+' p=' + str(correlation[1]) + \
                  '\n Linear regression coefficients: ' + 'm=' + str(m) + ' c=' +str(c)
        return correlation, [m, c]
