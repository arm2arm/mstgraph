#ifndef _METRIC_ND_
#define _METRIC_ND_

#include <boost/numeric/ublas/vector.hpp>
#include <vector>
#include "matrix_tools.h"

namespace Metrics{

	//////////////////////////////////////////////////////////////////////////////////
	//Compilcated Mahalanobis  distance for given particle from the group
	// let assume we have G group with NDxN components ahere ND is the number of the dimensions, 
	// and <G> is the means of the each component, so <G> is the ND vector.
	/** Usage Example:
	\code
	Metrics::MahalDistance<double> Mahal(uG1);
	ublas::vector<double> vec(2);
	vec(0)=uG1(0,0);vec(1)=uG1(1,0);
	cout<<"Pos:"<<vec<<endl;
	cout<<"Dist: "<<Mahal.doDist(vec)<<endl;
	\endcode
	**/
	template 
		<class T>	
	class MahalDistance {

		typedef boost::numeric::ublas::matrix<T> MatrixT;
		typedef boost::numeric::ublas::vector<T> VectorT;
	public:
		bool set_flag;
		MahalDistance():last_dist((T)0.0),set_flag(false){};		
		MahalDistance(MatrixT& G):last_dist((T)0.0){Init(G);};
		void Init(MatrixT& G){
			MatrixT Covar;
			Covarinace(G, Covar);//Get Covariance			
			InvertMatrix (Covar, InvCovar);
			for(size_t i=0;i<meanG.size();i++)
			meanV.push_back(meanG(i));
			detCovar=(T)determinant(Covar);
			//cout<<Covar<<endl;
			set_flag=true;
			}
		T doDist(float *V1,float *V2){
			boost::numeric::ublas::vector<T> p(meanG.size());
			for(size_t idx=0;idx<meanG.size();idx++)p(idx)=(V1[idx]-V2[idx]);
			last_dist=getDistance(p);
			return last_dist;
			};
		T doDist(float *V){
			boost::numeric::ublas::vector<T> p(meanG.size());
			for(size_t idx=0;idx<meanG.size();idx++)p(idx)=V[idx];
			last_dist=getDistance(p);
			return last_dist;
			};
		T doDist(std::vector<float> &V){
			boost::numeric::ublas::vector<T> p(meanG.size());
			for(size_t idx=0;idx<meanG.size();idx++)p(idx)=V[idx];
			last_dist=getDistance(p);
			return last_dist;
			};
		T doDist(VectorT &vin){
			last_dist=getDistance(vin);
			return last_dist;};
		inline	T getDet(){return detCovar;};
		T last_dist;
		T detCovar;
	protected:                  

		// this must be not directly accessible 
		// since we want to provide a rich set of distances	

		T getDistance(VectorT v1, VectorT v2)
			{
			return 1.0 ;
			};
		T getDistance(VectorT &vin)
			{
			
			vin-=meanG;
			VectorT ab=prod(vin,InvCovar );
			VectorT res=trans(vin);
			//cout<<ab<<endl;
			double dres=0.0;
			for(size_t i=0;i<res.size();i++)
				dres+=ab(i)*res(i);
			//cout<<InvCovar<<endl;
			//cout<<"Mean"<<meanG<<endl;
			return (T)dres;//
			};
		VectorT meanG;//this holds the means for each dimension
		MatrixT InvCovar;//this holds the means for each dimension
		std::vector<T> meanV;//just for debuging
		
		///////////////////////////// GET COVARIANCE of Matrix
		void Covarinace(const MatrixT& matrix, MatrixT& covMatrix)
			{
			const T ZERO(0);
			const T ONE(1);

			size_t maxRows = matrix.size2();
			size_t maxCols = matrix.size1();

			if (maxRows == 1) covMatrix = boost::numeric::ublas::zero_matrix<T>(1,1);
			else
				{

				GetMeanByColumns(matrix, meanG);
				MatrixT temp = matrix;
				for (size_t j=0; j<maxCols; ++j)
					{
					for (size_t i=0; i<maxRows; ++i) temp(j,i) -= meanG(j);
					}
				covMatrix = prod(temp,trans(temp)) * (ONE/T(maxRows-1));				
				}
			}

		void GetMeanByColumns(const MatrixT &m, VectorT & vout)
			{
			// loop over columns getting each row 
			vout.resize(m.size1());
			for (size_t i = 0; i < m.size1(); ++i) 
				{ 
				vout(i)=0.0;		
				for (size_t j = 0; j < m.size2(); ++ j)
					vout(i) += m(i,j);		
				vout(i)/=(T)m.size2();
				}
			}

		};	





	/////////////////////////////////////////////////////////////////////////////////

	}

#endif