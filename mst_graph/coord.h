#ifndef _MYCOORD_
#define _MYCOORD_


class CCoord;
typedef std::vector<CCoord> typeVecData;
class CCoord
	{
	public:
		CCoord(){};
		CCoord(int id_,MyFloat x_,MyFloat y_,MyFloat z_, MyFloat w_=0):
		id(id_), w(w_){
			pos[0]=x_;
			pos[1]=y_;
			pos[2]=z_;
			if(w==0.0)
				w=sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
			};
		MyFloat pos[3], vel[3];
		int id;
		MyFloat w;
		friend std::ostream& operator<<(std::ostream& os, const CCoord& v)
			{
			os.setf( std::ios::fixed, std::ios::floatfield ) ;
			os<<v.id;
			os.precision(2);
			os.width(8);
			os.setf( std::ios::fixed, std::ios::floatfield ) ;
			os<<" "<<v.pos[0];
			os.width(8);
			os<<"\t"<<v.pos[1];
			os.width(8);
			os<<"\t"<<v.pos[2];
			os.width(3);			
			os<<"\t w = "<<v.w<<std::endl;
			return os;
			}


	};


#endif