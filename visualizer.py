"""
visualizer.py:
Visualization tools for nsaba module.
Author: Torben Noto
"""
from nsaba import Nsaba
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D


class NsabaVisualizer(object):
    def __init__(self, nsaba_obj):
        self.no = nsaba_obj
    
    def visualize_ge(self, gene):
        for e in gene:
            if e in self.no.ge:
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                weights = self.no.ge[e]
                colors = cm.jet(weights/max(weights))
                color_map = cm.ScalarMappable(cmap=cm.jet)
                color_map.set_array(weights)
                fig.colorbar(color_map)

                x = self.no.aba['mni_coords'].data[:, 0]
                y = self.no.aba['mni_coords'].data[:, 1]
                z = self.no.aba['mni_coords'].data[:, 2]

                ax.scatter(x, y, z, c=colors, alpha=0.4)
            else:
                print 'Gene '+str(e) + ' has not been initialized. Use self.no.get_aba_ge([' + str(e) + '])'
        ax.set_title('Gene Expression of gene ID ' + str(gene))

        return fig

    def visualize_ns(self, term, points=200):
        if term in self.no.term:
            term_index = self.no.ns['features_df'].columns.get_loc(term)
            rand_point_inds = np.random.random_integers(0, len(np.squeeze(zip(self.no.ns['mni_coords'].data))), points)
            rand_points = np.squeeze(zip(self.no.ns['mni_coords'].data))[rand_point_inds]
            weights = []
            inds_of_real_points_with_no_fucking_missing_study_ids = []
            for rand_point in range(len(rand_points)):
                if len(self.no.coord_to_terms(rand_points[rand_point].astype(list))) > 0:
                    inds_of_real_points_with_no_fucking_missing_study_ids.append(rand_point_inds[rand_point])
                    weights.append(self.no.coord_to_terms(rand_points[rand_point].astype(list))[term_index])
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            colors = cm.jet(weights/max(weights))
            color_map = cm.ScalarMappable(cmap=cm.jet)
            color_map.set_array(weights)
            fig.colorbar(color_map)
            x = self.no.ns['mni_coords'].data[inds_of_real_points_with_no_fucking_missing_study_ids, 0]
            y = self.no.ns['mni_coords'].data[inds_of_real_points_with_no_fucking_missing_study_ids, 1]
            z = self.no.ns['mni_coords'].data[inds_of_real_points_with_no_fucking_missing_study_ids, 2]
        else:
            print 'Term '+term + ' has not been initialized. Use self.no.get_ns_act(' + term + ',thresh = 0.01)'
        ax.scatter(x, y, z, c=colors, alpha=0.4)
        ax.set_title('Estimation of ' + term)

        return fig

    def visualize_ns_ge(self, term, gene):
        for g in gene:
            if g in self.no.ge:
                if term in self.no.term:
                    ge_ns_mat = self.no.make_ge_ns_mat(term, gene)

                    # correlation
                    correlation = np.corrcoef(ge_ns_mat[:, 0], ge_ns_mat[:, 1])

                    # linear regression
                    X = np.vstack([ge_ns_mat[:, 0], np.ones(len(ge_ns_mat[:, 0]))]).T
                    m, c = np.linalg.lstsq(X, ge_ns_mat[:, 1])[0]

                    # print 'Correlation between ' + term + ' and gene number ' + str(gene)
                    # print correlation
                    # print 'Linear regression between ' + term + ' and gene number ' + str(gene) +' Slope =' + str(m) + ' y intercept = '+ str(c)
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    plt.plot(ge_ns_mat[:, 0], ge_ns_mat[:, 1], '.')
                    plt.plot([min(ge_ns_mat[:, 0]), max(ge_ns_mat[:, 0])], [m*min(ge_ns_mat[:, 0])+c, m*max(ge_ns_mat[:, 0])+c], 'r')
                    ax.set_xlabel(str(gene))
                    ax.set_ylabel(term)
                    return correlation, [m, c]
                else:
                    print 'Term '+term + ' has not been initialized. Use self.no.get_ns_act(' + term + ',thresh = 0.01)'
            else:
                print 'Gene '+str(g) + ' has not been initialized. Use self.no.get_aba_ge([' + str(g) + '])'

    def visualize_ns_ns(self, term1, term2):
        """Visualizing the relationship between two term vectors"""
        if term1 in self.no.term:
            if term2 in self.no.term:
                # correlation
                correlation = np.corrcoef(self.no.term[term1]['ns_act_vector'], self.no.term[term2]['ns_act_vector'])
                # linear regression
                regression_matrix = np.vstack([self.no.term[term1]['ns_act_vector'], np.ones(len(self.no.term[term1]['ns_act_vector']))]).T
                m, c = np.linalg.lstsq(regression_matrix, self.no.term[term2]['ns_act_vector'])[0]

                # plotting
                fig = plt.figure()
                ax = fig.add_subplot(111)
                plt.plot(self.no.term[term1]['ns_act_vector'], self.no.term[term2]['ns_act_vector'], '.')
                plt.plot([min(self.no.term[term1]['ns_act_vector']), max(self.no.term[term1]['ns_act_vector'])], [m*min(self.no.term[term2]['ns_act_vector'])+c, m*max(self.no.term[term2]['ns_act_vector'])+c], 'r')
                ax.set_xlabel(term1)
                ax.set_ylabel(term2)
                return correlation, [m, c]
            else:
                print 'Term '+term2 + ' has not been initialized. Use self.no.get_ns_act(' + term2 + ',thresh = 0.01)'
        else:
            print 'Term '+term1 + ' has not been initialized. Use self.no.get_ns_act(' + term1 + ',thresh = 0.01)'

    def visualize_ge_ge(self, genes):
        """Visualizing two gene vectors"""
        for gene in genes:
            if gene not in self.no.ge:
                print 'Gene '+str(g) + ' has not been initialized. Use self.no.get_aba_ge([' + str(g) + '])'
                return []
        # correlation
        correlation = np.corrcoef(self.no.ge[genes[0]], self.no.ge[genes[1]])
        # linear regression
        regression_matrix = np.vstack([self.no.ge[genes[0]], np.ones(len(self.no.ge[genes[0]]))]).T
        m, c = np.linalg.lstsq(regression_matrix, self.no.ge[genes[1]])[0]

        # plotting
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.plot(self.no.ge[genes[0]], self.no.ge[genes[1]], '.')
        plt.plot([min(self.no.ge[genes[0]]), max(self.no.ge[genes[0]])], [m*min(self.no.ge[genes[1]])+c, m*max(self.no.ge[genes[1]])+c], 'r')
        ax.set_xlabel(genes[0])
        ax.set_ylabel(genes[1])
        return correlation, [m, c]