#ifndef _METRIC_ND_
#define _METRIC_ND_
#undef BOOST_UBLAS_TYPE_CHECK
//we should define this to be able 
#define BOOST_UBLAS_TYPE_CHECK 0
#include <boost/numeric/ublas/vector.hpp>
#include <vector>
#include <algorithm>
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

		static bool set_flag;
		MahalDistance():last_dist((T)0.0),frac((T)0.5),NSmooth(5){	};		
		//! A constructor.
		/*!
		   It takes as a input the matrix of the D dimensional pointset.
		   Here the default values for frac and NSmooth
		   \sa frac  and NSmooth
		*/
		MahalDistance(MatrixT& G):last_dist((T)0.0),frac((T)0.5),NSmooth(5){Init(G);};
		//! A public member function.
		/*!
		   This used for setting the Covariance matrix smoothing.
		   \param frac the fraction of the particles to be used for smoothing
		   \param NSmooth number of iterations/shuffle to smooth covariance matrix  
		   \sa frac  and NSmooth
		*/
		
		void SetCovarianceSmoothing(T _frac,int _NSmooth){frac=_frac;NSmooth=_NSmooth;};
		void GetPartial(MatrixT &Gin, MatrixT &Gout, T frac)
			{
			
			std::random_shuffle(ind.begin(), ind.end());
			size_t n2=(size_t)(frac*Gin.size2())+1;
			size_t n1=Gin.size1();

			Gout.resize(Gin.size1(),n2);
			for(size_t i=0;i<n2;i++)
				for(size_t j=0;j<n1;j++)
				Gout(j,i)=Gin(j,ind[i]);

			}
		void Init(MatrixT& Gin){
			//cout<<endl;
			//cout<<G<<endl;

			MatrixT G;
			MatrixT Covar, partCovar;
			ind.resize(Gin.size2());
			for(size_t i=0;i<Gin.size2();i++)
				ind[i]=i;							

			GetPartial(Gin, G, frac);
			Covarinace(G, Covar);//Get Covariance			
			detCovar=(T)determinant(Covar);

			for(int i=0;i<NSmooth-1;i++)
				{
				GetPartial(Gin, G, frac);
				Covarinace(G, partCovar);//Get Covariance
				Covar+=partCovar;
				}
			Covar/=(T)NSmooth;
			InvertMatrix (Covar, InvCovar);
			
			detCovar=(T)determinant(Covar);
//			cout<<endl;
//			cout<<InvCovar<<endl;
			}
		T doDist(float *V1,float *V2){
			boost::numeric::ublas::vector<T> p(meanG.size());
			for(size_t idx=0;idx<meanG.size();idx++)p(idx)=(V1[idx]-V2[idx]);
			last_dist=getDistance(p);
			return last_dist;
			};
		T doEDist(float *V1,float *V2){
			last_dist=0.0;
			for(size_t idx=0;idx<3;idx++)
				last_dist+=(V1[idx]-V2[idx])*(V1[idx]-V2[idx]);
			last_dist=std::sqrt(last_dist);
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
		std::vector<int> ind;
	protected:                  
		T frac;// The fraction of neighbor  particles participating in to the Covariance Matrix smoothing
		int NSmooth;//	How many times select random fractions to smooth the covariance matrix?
		
		// this must be not directly accessible 
		// since we want to provide a rich set of distances	

		T getDistance(VectorT v1, VectorT v2)
			{
			return 1.0 ;
			};
		T getDistance(VectorT &vin)
			{
//			cout<<vin<<endl;
			//vin-=meanG;
			VectorT ab=prod(vin,InvCovar );
//			cout<<vin<<endl;
//			cout<<ab<<endl;
			double dres=0.0;
			for(size_t i=0;i<ab.size();i++)
				dres+=ab(i)*vin(i);
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
				MatrixT tr_temp=trans(temp);
				covMatrix = prod(temp,tr_temp) * (ONE/T(maxRows-1));				
				}
			}

		void GetMeanByColumns(const MatrixT &m, VectorT & vout)
			{
			// loop over columns getting each row 
			vout.resize(m.size1());
			meanV.clear();
			T mean=0.0;
			//cout<<m<<endl;
			for (size_t i = 0; i < m.size1(); ++i) 
				{ 
				mean=0.0;		
				for (size_t j = 0; j < m.size2(); ++ j)
					mean += m(i,j);		
				mean/=((T)m.size2());
				vout(i)=mean;
				meanV.push_back(mean);
			//	cout<<vout(i)<<endl;
				}
			}

		};	


template <class T>
bool MahalDistance<T>::set_flag=true;


	/////////////////////////////////////////////////////////////////////////////////

	}

#endif