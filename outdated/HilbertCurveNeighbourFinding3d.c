#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <sglib.h>

//#include <mpi.h>

#define min(a,b) a < b ? a : b
#define min4(a,b,c,d) min(min(a,b),min(c,d))
#define min8(a,b,c,d,e,f,g,h) min(min4(a,b,c,d),min4(e,f,g,h))
#define swap(a,b) typeof(a) c = a; a = b; b = c;

struct Point
{
	double x, y, z;
};
typedef struct Point Point;

Point makePoint(double a, double b,double c) {
	Point res;
	res.x = a;
	res.y = b;
	res.z = c;
	return res;
}

struct MDPoint {
	double *c;
};
typedef struct MDPoint MDPoint;

typedef unsigned int uint;

void MDPointScanf(MDPoint *p, uint dimensions) {
	int i;
	for(i=0;i<dimensions;i++) {
		scanf(" %lf",&p->c[i]);
	}
}

void MDPointInit(MDPoint *p,uint dimensions) {
	p->c = calloc(dimensions,sizeof(double));
}

MDPoint* MDPointCopy(MDPoint *other, uint dimensions) {
	MDPoint *new;
	MDPointInit(new,dimensions);
	memcpy(new->c,other->c,dimensions*sizeof(double));
	return new;
}

double MDPointHilbertPosition(MDPoint p, uint dimensions, double l) {
	if(l <= 0.0001) {
		return 0.0;
	}
	double suma = 0;
	double ile = 1;
	int i;
	for(i=dimensions-1;i>=0;i--) {
		printf("i=%d\n",i);
		if(p.c[i] < 0) {
			suma += l*ile;
			ile *= -1;
			p.c[i] *= (-1);
		}
		l/=2;
	}
	printf("l:%f, suma :: %f\n",l,suma);
	for(i=0;i<dimensions;i++) {
		printf("p.c[%d] == %f\n",i,p.c[i]);
		p.c[i] -= l;
	}
	return suma + MDPointHilbertPosition(p,dimensions,l)*ile;
}

struct Pair {
	void *first;
	void *second;
};
typedef struct Pair Pair;

#define AllocatePair(obj,x,y) \
obj->first = malloc(sizeof(x)); \
obj->second = malloc(sizeof(y));

#define SetPair(obj,x,y) \
AllocatePair(obj,x,y); \
*(typeof(x)*)obj->first = x; \
*(typeof(y)*)obj->second = y; 

#define MakePair()\
malloc(sizeof(Pair))

bool pointComparator(Point a, Point b) {
	if(a.x != b.x)
		return a.x < b.x;
	if(a.y != b.y)
		return a.y < b.y;
	return a.z < b.z;
}

Point mirror(int a, int b, int c, Point p) {return makePoint(a * p.x, b * p.y, c * p.z); }

int sgn(double x)
{
        if(x < 0) return -1;
        if(x >= 0) return 1;
}

Point trans(Point p)
{
        p.x -= 0.5 * sgn(p.x);
        p.y -= 0.5 * sgn(p.y);
        p.z -= 0.5 * sgn(p.z);
        p.x *= 2.0; p.y *= 2.0; p.z *= 2.0;
        return p;
}

Point rotate(Point p, char a, int k)
{
        if(a == 'z')
        {
                swap(p.x, p.y);
                p.y *= k;
                p.x *= -k;
        }
        else if(a == 'x')
        {
                swap(p.y, p.z);
                p.y *= -k;
                p.z *= k;
        }
        else if(a == 'y')
        {
                swap(p.x, p.z);
                p.x *= k;
                p.z *= -k;
        }
        return p;
}



#define EPS 1e-12	

inline double HilbertPosNonRecursive(Point p) {
	double l = 1.0;
	double sum = 0;
	double znak = 1;
	while( l >= EPS) {
		if(p.z > 0) {
			sum += l*znak;
			znak *= -1;
			p.z = -p.z;
		} else {
			l/=8.0;
		        if(p.x >= 0 && p.y >= 0) p = rotate( rotate(trans(p), 'z' , 1), 'y', -1);
		        else if(p.x >= 0 && p.y < 0) { 
				sum += l * znak;
				p = rotate( rotate(trans(p), 'y' , 1), 'z', -1);
			}
        		else if(p.x < 0 && p.y < 0)  {
				sum += l*2.0 * znak;
				p = rotate( rotate(trans(p), 'y' , 1), 'z', -1);
			}
		        else if(p.x < 0 && p.y >= 0) {
				sum += l*3.0 * znak;
				p = rotate( rotate(trans(p), 'z' , 1), 'z', 1); 
			}
		}
	}	
	return sum;
}

#define HilbertPos(p) RealHilbertPos(p,1.0)
double RealHilbertPos(Point p, double l)
{
        if(l < EPS) 
		return 0.0;
        if(p.z > 0) 
		return l - RealHilbertPos(makePoint(p.x, p.y, -p.z), l);
        if(p.x >= 0 && p.y >= 0) 
		return RealHilbertPos( rotate( rotate(trans(p), 'z' , 1), 'y', -1), l / 8.0 );
        if(p.x >= 0 && p.y < 0) 
		return l / 8.0 + RealHilbertPos( rotate( rotate(trans(p), 'y' , 1), 'z', -1), l / 8.0 );
        if(p.x < 0 && p.y < 0) 
		return l / 4.0 + RealHilbertPos( rotate( rotate(trans(p), 'y' , 1), 'z', -1), l / 8.0 );
        if(p.x < 0 && p.y >= 0) 
		return 3.0 * l / 8.0 + RealHilbertPos( rotate( rotate(trans(p), 'z' , 1), 'z', 1), l / 8.0 );
}


inline bool intersects(Point u, double d1, Point v, double d2) // if square which center is in u intersects sqare which center is in v (d1 and d2 are halfs of square sides)
{
	if(abs(u.x - v.x) <= d1 + d2 && abs(u.y - v.y) <= d1 + d2 && abs(u.z - v.z) <= d1 + d2) return true;
	return false; 
}

inline bool contains(Point u, double d1, Point v, double d2) // if square which center is in u contains square which center is in v (d1 and d2 are halfs of square sides)
{
	return 
	u.x - d1 >= v.x - d2 && 
	u.x + d1 <= v.x + d2 && 
	u.y - d1 >= v.y - d2 && 
	u.y + d1 <= v.y + d2 &&
	u.z - d1 >= v.z - d2 && 
	u.z + d1 <= v.z + d2;
}

void query(Point u, double a, Point p, double r, double l, Pair *res, int *ressize) // query adds intervals on Hilbert Curve to res, where neighbours of p may appear
{
	if(!intersects(u, a, p, r)) return;
	if(contains(u, a, p, r) || l < 0.01) 
	{
		double c = min8(
		 HilbertPosNonRecursive(makePoint(u.x + a / 2.0, u.y + a / 2.0, u.z + a/2.0)),
		 HilbertPosNonRecursive(makePoint(u.x - a / 2.0, u.y + a / 2.0, u.z + a/2.0)), 
		 HilbertPosNonRecursive(makePoint(u.x - a / 2.0, u.y - a / 2.0, u.z + a/2.0)), 
		 HilbertPosNonRecursive(makePoint(u.x + a / 2.0, u.y - a / 2.0, u.z + a/2.0)),
		 HilbertPosNonRecursive(makePoint(u.x + a / 2.0, u.y + a / 2.0, u.z - a/2.0)),
		 HilbertPosNonRecursive(makePoint(u.x - a / 2.0, u.y + a / 2.0, u.z - a/2.0)), 
		 HilbertPosNonRecursive(makePoint(u.x - a / 2.0, u.y - a / 2.0, u.z - a/2.0)), 
		 HilbertPosNonRecursive(makePoint(u.x + a / 2.0, u.y - a / 2.0, u.z - a/2.0))
		); 
		//res[(*ressize)++] = makePair(c,c+l)
	} else {
		query(makePoint(u.x + a / 2.0, u.y + a / 2.0, u.z + a/2.0),a/2.0, p, r,  l/8.0, res, ressize);
		query(makePoint(u.x - a / 2.0, u.y + a / 2.0, u.z + a/2.0),a/2.0, p, r,  l/8.0, res, ressize);
		query(makePoint(u.x + a / 2.0, u.y + a / 2.0, u.z - a/2.0),a/2.0, p, r,  l/8.0, res, ressize);
		query(makePoint(u.x - a / 2.0, u.y + a / 2.0, u.z - a/2.0),a/2.0, p, r,  l/8.0, res, ressize);
		query(makePoint(u.x + a / 2.0, u.y - a / 2.0, u.z + a/2.0),a/2.0, p, r,  l/8.0, res, ressize);
		query(makePoint(u.x - a / 2.0, u.y - a / 2.0, u.z + a/2.0),a/2.0, p, r,  l/8.0, res, ressize);
		query(makePoint(u.x + a / 2.0, u.y - a / 2.0, u.z - a/2.0),a/2.0, p, r,  l/8.0, res, ressize);
		query(makePoint(u.x - a / 2.0, u.y - a / 2.0, u.z - a/2.0),a/2.0, p, r,  l/8.0, res, ressize);
	}
}



struct Assign {
	double p;
	double k;
	int process;
};

typedef struct Assign Assign;

Assign makeAssign(double a, double b, int c) {
	Assign r;
	r.p = a;
	r.k = b;
	r.process = c;
	return r;
}

int assignComparator(Assign a, Assign b) {
	if(a.p != b.p)
		return a.p - b.p;
	else
		return a.k - b.k;
}

/*struct Mylist {
	Pair p;
	Mylist* next_ptr;
};

Assign *assignements_list;
Mylist *worek;

void dodaj(Pair *res, int ressize) {
	int i;
	for(i=0;i<ressize;i++) {
		Mylist *l = malloc(sizeof(Mylist));
		SGLIB_LIST_ADD(struct Mylist,worek,l,next_ptr);		
	}		
}



Assign* cores;

void dodaj_przedzial(double p, double k, int coresNum) {
	vector <assignement>::iterator it = lower_bound(cores.begin(), cores.end(), assignement(p, 0, 0));
	int left,right,middle;
	for(left=0,right=0;left<right;) {
		middle = (left + right)/2;
		if(cores[middle].p >= p)
			right = middle;
		else
			left = middle+1;
	}
	
}

int main () {return 0;}



const double EPS420 = 1e-9;

bool blisko(double a, double b) {
	if(a > b)
		swap(a,b);
	if(a + EPS420 > b)
		return true;
	return false;
}
*/

int main(int argc, char *argv[]) {
	MDPoint p;
	MDPointInit(&p,3);
	MDPointScanf(&p,3);
	printf("%f",MDPointHilbertPosition(p,3,1.0));
	return 0;		
}
	/*MPI_Init(&argc, &argv);
	vector< pivpdd > odpowiedzi;
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	double radius; // radius - specify how far a point must be to be a neighbour
	int pointCount; // number of points assigned to the process

	vector<Point> coordinates;
	vector<double> curvePos;
	cores = calloc(size,sizeof(Assign));
	if(rank == 0) {
		FILE *fp;
		fp = fopen("tests/HilbertInput3D");
		int n; cin >> n;
		vector <Point> points(n);
		vector <pdi> HilPos(n);
		for(int i=0; i < n; i++)
		{
			cin >> points[i];
			HilPos[i] = make_pair( HilbertPosNonRecursive(points[i]), i );
		}
		
		sort(HilPos.begin(), HilPos.end());
		typedef vector<int> vi;
		vector< vi >  procki(size);
		double *intervals = new double [2*size+1];
		for(int i=0;i<size;i++) {
			intervals[2*i] = 2;
			intervals[2*i+1] = -1;
		}
		int czesc = n/size;
		int ile = 0;
		
		for(int i = 0; i < n; ++i)  {
			int akt = i*size/n;
			procki[akt].push_back(HilPos[i].second);
			//cout << "do procka "  << akt << " przyporzodkowuje " << points[HilPos[i].second].x << ", "  << points[HilPos[i].second].y << endl;
			intervals[2*akt] = min(intervals[2*akt], HilPos[i].first);
			intervals[2*akt+1] = max(intervals[2*akt+1], HilPos[i].first);
		}
		intervals[2*size] = radius;
	 	MPI_Bcast( intervals, 2*size+1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

		int *arr = new int[size];
		for(int i=0;i<size;i++) {
			arr[i] = procki[i].size();
			cout << "procesor " << i << " ma " << arr[i] << endl;
		}
		pointCount = arr[0];
		//cout << endl;
		int *recvBuffer = new int[1];
		MPI_Scatter(arr, 1, MPI_INT, recvBuffer, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		double *wspolrzedne = new double [n*4];
		for(int i=0;i<n*4;i+=4) {
			wspolrzedne[i] = points[HilPos[i/4].second].x;
			wspolrzedne[i+1] = points[HilPos[i/4].second].y;
			wspolrzedne[i+2] = points[HilPos[i/4].second].z;
			wspolrzedne[i+3] = HilPos[i/4].first;
		}
		
		int catfish = arr[0]*4;
		for(int i=1;i<size;i++) {
			MPI_Send(wspolrzedne+catfish, arr[i]*4, MPI_DOUBLE, i,1,MPI_COMM_WORLD);
			catfish += arr[i]*4;
		}
		for(int i=0;i<arr[0];i++) {
			coordinates.push_back(Point(wspolrzedne[i*4],wspolrzedne[i*4+1], wspolrzedne[i*4+2]));
			curvePos.push_back(wspolrzedne[i*4+3]);
		}
		for(int i=0;i<size;i++) {
			cores.push_back(assignement(intervals[i*2],intervals[i*2+1],i));
		}

		MPI_Barrier(MPI_COMM_WORLD);
		//cout << "kończy procesor " << rank << endl;
	}
	else
	{	
		double *inter = new double[2*size+1];
	 	MPI_Bcast( inter, 2*size+1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
		radius = inter[2*size];
		int* recvBuffer = new int[1];
		MPI_Scatter(NULL,1,MPI_INT,recvBuffer,1,MPI_INT,0,MPI_COMM_WORLD);
		//cout << "Jestem procesorem #" << rank << " i dostalem " << recvBuffer[0] << endl;
		double *coor = new double[recvBuffer[0]*4];
		pointCount = recvBuffer[0];
		MPI_Status stat;
		MPI_Recv(coor, recvBuffer[0]*4, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &stat);
		for(int i=0;i<size;i++) {
			cores.push_back(assignement(inter[i*2],inter[i*2+1],i));
		}
		for(int i=0;i<recvBuffer[0];i++) {
			coordinates.push_back(Point(coor[i*4], coor[i*4+1], coor[i*4+2]));
			curvePos.push_back(coor[i*4+3]);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		//cout << "kończy procesor " << rank << endl;
	}
	vector <pdd> res;
	for(int i=0;i<pointCount;i++) {
		res.clear();
		query(Point(0.0, 0.0, 0.0), 1.0, coordinates[i], radius, 1.0, res);
		odpowiedzi.push_back(pivpdd(i,res));	
		if(res.size() == 0) {
			cout << "jest zle w " << coordinates[i] << endl;
		}
		assert(res.size()>0);
	}
	
	assert(is_sorted(cores.begin(),cores.end()));
	sort(cores.begin(), cores.end());
	
	
	for(int i=0;i<(int)odpowiedzi.size();i++) {
		dodaj(odpowiedzi[i].second);	
	}

	assert(worek.size() > 0);
	
	vector<pdi> events;
	for(int i=0;i<(int)worek.size();i++) {
		events.push_back(pdi(worek[i].first,1));
		events.push_back(pdi(worek[i].second+1e-12,-1));			
	}
	sort(events.begin(),events.end());
	
	vector<pdd> wspolne;
	int suma = 1;
	double last = events[0].first;
	
	for(int i=1;i<(int)events.size();i++) {
			if(suma > 0) {
				wspolne.push_back(pdd(last,events[i].first-1e-12));
				while(wspolne.size() > 1 && blisko(wspolne[wspolne.size()-2].second,wspolne[wspolne.size()-1].first)) {
					wspolne[wspolne.size()-2].second = wspolne[wspolne.size()-1].second;
					wspolne.pop_back();
				}
			}
			suma += events[i].second;
			last = events[i].first;	
	}
	
	for(int i=0;i<(int)wspolne.size();i++) {
		assert(wspolne[i].first-1e-12 <= wspolne[i].second);
		dodaj_przedzial(wspolne[i].first, wspolne[i].second);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	int* ile_od_sasiada = new int[size];
	for(int i=0;i<size;i++)
		ile_od_sasiada[i] = 0;
	for(int i=0;i<(int)assignements.size();i++) {
		if(assignements[i].procesor != rank)
			ile_od_sasiada[assignements[i].procesor]++;
	}
	int* ile_do_mnie = new int[size];
	MPI_Alltoall(ile_od_sasiada, 1, MPI_INT,ile_do_mnie, 1, MPI_INT,MPI_COMM_WORLD);

	
	int ptr = 0;
	MPI_Request *sendingRequest = new MPI_Request[size];
	for(int i=0;i<size;i++) if(i!=rank && ile_od_sasiada[i] != 0) {
		double *intervalsPls = new double [ile_od_sasiada[i]*2];
		int ilePls = 0;
		while(ptr < assignements.size() && assignements[ptr].procesor <= i) {
			if(assignements[ptr].procesor == i) {
				intervalsPls[ilePls] = assignements[ptr].p;
				intervalsPls[ilePls+1] = assignements[ptr].k;
				ilePls+=2;
			}
			ptr++;
		}
		MPI_Isend(intervalsPls, 2*ile_od_sasiada[i], MPI_DOUBLE,i,2,MPI_COMM_WORLD,&sendingRequest[i]);
	}


	double **intervalsSir = new double* [size];
	MPI_Request *mpir2 = new MPI_Request[size];
	for(int i=0;i<size;i++) {
		if(i==rank or ile_do_mnie[i] == 0) {
			intervalsSir[i] = NULL;
			continue;
		}
		intervalsSir[i] = new double[ile_do_mnie[i]*2];
		MPI_Irecv(intervalsSir[i], ile_do_mnie[i]*2, MPI_DOUBLE, i,2, MPI_COMM_WORLD, &mpir2[i]);
	}
	MPI_Status *stat = new MPI_Status[size*2];
	for(int i=0;i<size;i++) if(i!=rank && ile_od_sasiada[i] != 0)
		MPI_Wait(&sendingRequest[i],&stat[i]);
	for(int i=0;i<size;i++) if(i!=rank && ile_do_mnie[i] != 0)
		MPI_Wait(&mpir2[i],&stat[size+i]);
	MPI_Barrier(MPI_COMM_WORLD);


	int *ile_punktow_ma_dla_mnie_sasiad = new int[size];
	int *ilejamamdla = new int[size];
	
	vector<Point> *wyn = new vector<Point>[size];
	vector<double> *wyn77 = new vector<double>[size];
	for(int i = 0; i < size; ++i) if(intervalsSir[i] != NULL) {
		wyn[i].clear();
		wyn77[i].clear();
		for(int j = 0; j < 2*ile_do_mnie[i]; j+= 2) {
			pdd f = pdd(intervalsSir[i][j],intervalsSir[i][j+1]);
			vector<double>::iterator it = lower_bound(curvePos.begin(),curvePos.end(),f.first);
			while(it != curvePos.end() && *it <= f.second) {
				int kt = it-curvePos.begin();
				wyn[i].push_back(coordinates[kt]);
				wyn77[i].push_back(*it);
				it++;
			}
		}
	}

	for(int i=0;i<size;i++)
		ilejamamdla[i] = wyn[i].size();
	MPI_Alltoall(ilejamamdla, 1, MPI_INT,ile_punktow_ma_dla_mnie_sasiad, 1, MPI_INT,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	

	MPI_Request *neighbourSendMPIRequest = new MPI_Request[size];

	double **buf = new double*[size];
	for(int i=0;i<size;i++) if(i!=rank && wyn[i].size() > 0) {
		buf[i] = new double[wyn[i].size()*4];
		for(int k=0;k<(int)wyn[i].size();k++) {
			buf[i][k*4] = wyn[i][k].x;
			buf[i][k*4+1] = wyn[i][k].y;
			buf[i][k*4+2] = wyn[i][k].z;
			buf[i][k*4+3] = wyn77[i][k];
		}
		MPI_Isend(buf[i], wyn[i].size()*4, MPI_DOUBLE, i,7,MPI_COMM_WORLD, &neighbourSendMPIRequest[i]);
	}
	
	MPI_Request *neighbourRecvMPIRequest = new MPI_Request[size];

	vector<pdp> box;
	
	
	double **buff = new double*[size];

	for(int i=0;i<size;i++) if(i!=rank && ile_punktow_ma_dla_mnie_sasiad[i] > 0) {
		buff[i] = new double[ile_punktow_ma_dla_mnie_sasiad[i]*4];
		MPI_Irecv(buff[i], ile_punktow_ma_dla_mnie_sasiad[i]*4, MPI_DOUBLE,i,
		              7, MPI_COMM_WORLD, &neighbourRecvMPIRequest[i]);
	}
	


	for(int i=0;i<size;i++) if(i!=rank && wyn[i].size() > 0)  {
		MPI_Wait(&neighbourSendMPIRequest[i],&stat[i]);
	}
	for(int i=0;i<size;i++) if(i!=rank && ile_punktow_ma_dla_mnie_sasiad[i] > 0) {
		MPI_Wait(&neighbourRecvMPIRequest[i],&stat[i+size]);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	

	for(int i=0;i<size;i++) {
		for(int j=0;j<ile_punktow_ma_dla_mnie_sasiad[i]*4;j+=4) {
			box.push_back(pdp(buff[i][j+3],Point(buff[i][j],buff[i][j+1], buff[i][j+2])));
		}
	}

	sort(box.begin(),box.end());

	for(int i=0;i<box.size();i++) {
		cout << "jestem " << rank << " Mam sasiada!  : " << box[i].second << endl;
	}

	
	//cout << "Jestem procesem " << rank << " i kończę\n";
	MPI_Finalize();
	return 0;

}*/ 

