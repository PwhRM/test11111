#include <iostream>
//for ploting
//gnuplot include must be declared in the first
//#include "gnuplot_i.hpp"

#include <armadillo>
//cluster data simulation
#include "clust_sim.hpp"
//make Euclidean distance matrix
#include "dist_mat_maker.hpp"

/*
//hist_eq just for testing
#include "hist_eq.hpp"

//inImDynM
#include "inImDynM.hpp"

//DBSCAN expand cluster
#include "ExpandCluster.hpp"
 */

//DSets-DBSCAN
#include "DSets_DBSCAN.hpp"

//diplay image using opencv
//#include <opencv2/core/core.hpp>
//#include <opencv2/highgui/highgui.hpp>



using namespace std;
using namespace arma;

//using std::cout;
//using std::endl;

int main(int argc, const char **argv) {
    //generate random cluster data
    //simulation parameters
    /*
    int n_clust = 1, n_point = 10, n_noise = 0;
    double r = 25, height = 1000, width = 1000;
    mat clust_data;
    clust_data = clust_sim(n_clust, n_point, n_noise, r, height, width);
    */
    mat clust_data = {{1595, 1671, 1623, 1580, 1694, 1648, 1677, 1658, 1688, 1598},
        {217, 289, 316, 233, 234, 247, 322, 248, 247, 262}};
    clust_data = clust_data.t();
    //distance matrix of cluster data
    mat dist_mat = dist_mat_maker(clust_data.cols(0, 1));
    
    vec clust(dist_mat.n_cols);
    
    clust = DSets_DBSCAN(dist_mat);
    /*
    double sigma = 10 * mean(mean(dist_mat));
    //test hist_eq function
    mat A = exp(- dist_mat / sigma);
    A.diag() = A.diag() - A.diag();
    int hist_size = 50;
    //histogram equalization
    A = hist_eq(A, hist_size);
    //A.diag() = A.diag() - A.diag();
    
    
    //test other functio in DSets-DBSCAN
    //number of points in the data set
    int n_points = int(A.n_rows);
    //index of points in the orignal data set
    //0 value will be ocupied as visited points so the index start from 1
    vec d(A.n_rows);
    for (int i = 0; i < n_points; i++) {
        d(i) = i + 1;
    }
    
    //initialize cluster ID and other parameters
    int cid = 0;
    int minPts = 4;
    int maxIters = 1000;
    double supportThreshold = 1e-4;
    double precision = 1e-8;
    //n dimensional simplex as x0
    mat x0 = ones(n_points, 1);
    x0 = x0 / sum(sum(x0));
    //container for cluster ID
    vec clust = zeros(n_points) - 1;
    //while loop control
    int count = n_points;
    int iter_count = 0;
    while (any(d) && count > 0) {
        cid++;
        //create temporal d, A and dist_mat to be used in the loop
        //d(i) = 0: point i is visited
        vec d_tmp = d(find(d > 0));
        cout << "d_tmp: " << d_tmp.t() << endl;
        cout << "clust: " << clust.t() << endl;
        //number of un-visited points
        int n_points_tmp = int(d_tmp.n_elem);
        //debug
        count = n_points_tmp;
        iter_count++;
        //cout << d_tmp.t() << endl;
        cout << "number of iteration: " << iter_count << ", count = " << count << endl;
        
        mat A_tmp(n_points_tmp, n_points_tmp),
        dist_mat_tmp(n_points_tmp, n_points_tmp);
        for (int i = 0; i < n_points_tmp; i++) {
            for (int j = 0; j < n_points_tmp; j++) {
                A_tmp(i, j) = A(d_tmp(i) - 1, d_tmp(j) - 1);
                dist_mat_tmp(i, j) = dist_mat(d_tmp(i) - 1, d_tmp(j) - 1);
            }
        }
        //get current x from current d
        //the size of d is changing, so size of x is changing too
        mat x(d_tmp.n_elem, 1);
        x = ones(d_tmp.n_elem, 1);
        x = x / arma::sum(arma::sum(x));
        //extract dominant from DSets method
        x = inImDynM(A_tmp, x, precision, maxIters);
        //note id_A will be vec instead of matrix
        //x_good_ind: index of good x in vector x
        uvec x_good_ind = find(x > supportThreshold);
        uvec id_A = x_good_ind;
        //assign clust id
        //id_clust - 1 is the actual index in clust
        vec id_clust = d_tmp.elem(x_good_ind);
        for (int i = 0; i < int(id_clust.n_elem); i++) {
            clust(id_clust(i) - 1) = cid;
        }
        //use DBSCAN to expand dominant
        double Eps = 0;
        for (int i = 0; i < int(id_A.n_elem); i++) {
            //dist_tmp should be the distance map in dist_mat, not dist_mat_tmp
            vec dist_tmp = dist_mat_tmp.col(id_A(i));
            dist_tmp = arma::sort(dist_tmp);
            Eps = max(Eps, dist_tmp(std::min(minPts, int(dist_tmp.n_elem)) - 1));
        }
        uvec visit_list = x_good_ind;
        //debug
        cout << "cluster index by DSets: " << endl;
        cout << id_clust.t() - 1 << endl;
        for (int i = 0; i < int(x_good_ind.n_elem); i++) {
            uword pt = visit_list(i);
            //pt will be pt - 1 inside the ExpandCluster function
            clust = ExpandCluster(dist_mat_tmp, pt, cid, Eps, minPts, clust, d_tmp);
        }
        //get new xid
        vec d_clust = d(find(clust == cid));
        cout << "cluster index after DBSCAN: " << endl;
        cout << d_clust.t() - 1 << endl;
        //assigan visited data points in d as 0
        for (int i = 0; i < d_clust.n_elem; i++) {
            d(d_clust(i) - 1) = 0;
        }
    }
     */
    cout << clust << endl;
    
    
    
    
    /*
    //plot data for tets
    Gnuplot::set_terminal_std("qt");
    
    Gnuplot g1("lines");
    cout << "*** plotting slopes" << endl;
    g1.set_title("old cdf");
    //convert arma mat to std::vector
    vector<double> x1 = conv_to<vector<double>>::from(x.col(0));
    vector<double> y1 = conv_to<vector<double>>::from(cdf.col(0));
    g1.plot_xy(x1, y1);
    g1.set_legend();
    */
    /*
    //display dist_mat and similarity matrix for test
    cv::Mat image_test(int(A.n_rows), int(A.n_cols), CV_64FC1, A.memptr());
    cv::namedWindow("display", cv::WINDOW_AUTOSIZE);
    cv::imshow("display", image_test);
    */
    
    //using input as waitkey()
    int waittest;
    cin >> waittest;
    return 0;
}
