
#include <bits/stdc++.h>
#include <mpi.h>
using namespace std;


struct Point // Point structure
{
	double x, y;
	Point() {}
	Point(double a, double b) {x = a; y = b;}	
	void rotate(int k) {
		swap(x, y);
		y *= -k;
		x *= k;
	}
	Point operator+(Point p) {return Point(x + p.x, y + p.y); } // addition operator 
	Point operator*(double k) {return Point(k * x, k * y); } // multiplication operator
	double distance(Point p) {return sqrt((x-p.x)*(x-p.x) + (y-p.y)*(y-p.y));} // distance counting in Euclidean 
};

bool operator < (Point a, Point b) { // nothing relevant - just for sorting purposes
	return a.x < b.x;
}

Point rotate(Point p, int k) { p.rotate(k); return p; }
Point mirrorX(Point p) {return Point(-p.x, p.y); }

#define EPS 1e-4	

inline double HilbertPos(Point p, double l = 1.0) // Hilbert Curve Position Finding, returns value in [0,1] for given point
{
	if(l < EPS) return 0.0;
	if(p.x <= 0 && p.y <= 0) {return l / 4.0 + HilbertPos((p + Point(0.5, 0.5)) * 2.0, l / 4.0); }
	if(p.x <= 0 && p.y > 0) {return HilbertPos(mirrorX(rotate((p + Point(0.5, -0.5)) * 2.0, 1)), l / 4.0); }
	if(p.x > 0 && p.y <= 0) {return l / 2.0 + HilbertPos((p + Point(-0.5, 0.5)) * 2.0, l / 4.0); }
	if(p.x > 0 && p.y > 0) {return 3.0 * l / 4.0 + HilbertPos(mirrorX(rotate((p + Point(-0.5, -0.5)) * 2.0, -1)), l / 4.0); }
	return 0.0;
}

// some useful typedefs
typedef pair <int, Point> PIP;
typedef pair <int,int> pii;
typedef pair<double, Point> pdp;
typedef pair<double,double> pdd;
typedef pair< int, vector<pdd> > pivpdd;
typedef pair< double, int> pdi;

inline bool intersects(Point u, double d1, Point v, double d2) // if square which center is in u intersects sqare which center is in v (d1 and d2 are halfs of square sides)
{
	if(abs(u.x - v.x) <= d1 + d2 && abs(u.y - v.y) <= d1 + d2) return true;
	return false; 
}

inline bool contains(Point u, double d1, Point v, double d2) // if square which center is in u contains square which center is in v (d1 and d2 are halfs of square sides)
{
	return u.x - d1 >= v.x - d2 && u.x + d1 <= v.x + d2 && u.y - d1 >= v.y - d2 && u.y + d1 <= v.y + d2;
}

inline double min(double a, double b, double c, double d) { // min function for 4 numbers
	return min(min(a,b),min(c,d));
}

void query(Point u, double a, Point p, double r, double l, vector <pdd> &res) // query adds intervals on Hilbert Curve to res, where neighbours of p may appear
{
	if(!intersects(u, a, p, r)) return;
	if(contains(u, a, p, r) or l < 0.01) 
	{
		double c = min(
		 HilbertPos(Point(u.x + a / 2.0, u.y + a / 2.0)),
		 HilbertPos(Point(u.x - a / 2.0, u.y + a / 2.0)), 
		 HilbertPos(Point(u.x - a / 2.0, u.y - a / 2.0)), 
		 HilbertPos(Point(u.x + a / 2.0, u.y - a / 2.0))
		); 
		res.push_back(pdd(c, c + l));
	} else {
		query(Point(u.x + a / 2.0, u.y + a / 2.0), a / 2.0, p, r, l / 4.0, res);
		query(Point(u.x - a / 2.0, u.y + a / 2.0), a / 2.0, p, r, l / 4.0, res);
		query(Point(u.x - a / 2.0, u.y - a / 2.0), a / 2.0, p, r, l / 4.0, res);
		query(Point(u.x + a / 2.0, u.y - a / 2.0), a / 2.0, p, r, l / 4.0, res);
	}
}



struct assignement {
	double p;
	double k;
	int procesor;
	assignement () {}
	assignement (double a, double b, int c) {p = a, k = b, procesor = c;}	
};

bool operator < (assignement x, assignement y) {
	return pdd(x.p,x.k) < pdd(y.p,y.k);
}

vector<assignement> assignements;

vector<pdd> worek;

inline void dodaj( vector<pdd> &res) {
	for(int i=0;i<(int)res.size();i++) {
		worek.push_back(res[i]);
	}		
}

/*double force_van_der_vaals(Point x, Point y) // van der vaals physics
{
	double r = x.distance(y);
	double a = 1.0;
	double b = 1.0;
	return a / pow(r, 13.0) - b / pow(r, 7.0);
}*/

vector<assignement> cores;

void dodaj_przedzial(double p, double k) {
	vector <assignement>::iterator it = lower_bound(cores.begin(), cores.end(), assignement(p, 0, 0));
	
	if(it != cores.begin()) {
		it--;
		if(p < it -> k) assignements.push_back(assignement(p, min(it -> k, k), it -> procesor));
		it++;
	}
	
	while (it != cores.end() && it -> p <= k) {
		assignements.push_back(assignement(it -> p, min(it -> k, k), it -> procesor));
		it++;
	}		
}

const double EPS420 = 1e-9;

bool blisko(double a, double b) {
	if(a > b)
		swap(a,b);
	if(a + EPS420 > b)
		return true;
	return false;
}

double radius; // radius - specify how far a point must be to be a neighbour
int pointCount; // number of points in single process

int main(int argc, char *argv[]) {
	vector< pivpdd > odpowiedzi;
	freopen("HilbertInput","r",stdin);
	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	vector<Point> coordinates;
	vector<double> curvePos;

	if(rank == 0) {
		cin >> radius;
		int n; cin >> n;
		vector <Point> points(n);
		vector <pdi> HilPos(n);
		for(int i=0; i < n; i++)
		{
			cin >> points[i].x >> points[i].y;
			HilPos[i] = make_pair( HilbertPos(points[i], 1.0), i );
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
		/*for(int i = 0; i < 2 * size; ++i) cout << intervals[i] << ' ';
		cout << endl;*/

		int *arr = new int[size];
		for(int i=0;i<size;i++) {
			arr[i] = procki[i].size();
			//cout << "procesor " << i << " ma " << arr[i] << endl;
		}
		pointCount = arr[0];
		//cout << endl;
		int *recvBuffer = new int[1];
		MPI_Scatter(arr, 1, MPI_INT, recvBuffer, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		double *wspolrzedne = new double [n*3];
		for(int i=0;i<n*3;i+=3) {
			wspolrzedne[i] = points[HilPos[i/3].second].x;
			wspolrzedne[i+1] = points[HilPos[i/3].second].y;
			wspolrzedne[i+2] = HilPos[i/3].first;
		}
		
		int catfish = arr[0]*3;
		for(int i=1;i<size;i++) {
			MPI_Send(wspolrzedne+catfish, arr[i]*3, MPI_DOUBLE, i,1,MPI_COMM_WORLD);
			catfish += arr[i]*3;
		}
		for(int i=0;i<arr[0];i++) {
			coordinates.push_back(Point(wspolrzedne[i*3],wspolrzedne[i*3+1]));
			curvePos.push_back(wspolrzedne[i*3+2]);
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
		/*for(int i=0;i<size*2;i+=2) {
			cout << "(" << inter[i] << ", " << inter[i+1] << ") ";
		}
		cout << endl;*/
		radius = inter[2*size];
		int* recvBuffer = new int[1];
		MPI_Scatter(NULL,1,MPI_INT,recvBuffer,1,MPI_INT,0,MPI_COMM_WORLD);
		//cout << "Jestem procesorem #" << rank << " i dostalem " << recvBuffer[0] << endl;
		double *coor = new double[recvBuffer[0]*3];
		pointCount = recvBuffer[0];
		MPI_Status stat;
		MPI_Recv(coor, recvBuffer[0]*3, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &stat);
		//cout << "Jestem procesorem #" << rank << " i dostalem wspolrzedne :\n";
		/*for(int i=0;i<recvBuffer[0]*3;i+=3) {
			cout << "Jestem " << rank <<  "(" << coor[i] << ", " << coor[i+1] << ") -> " << coor[i+2] << endl;
		}*/
		for(int i=0;i<size;i++) {
			cores.push_back(assignement(inter[i*2],inter[i*2+1],i));
		}
		for(int i=0;i<recvBuffer[0];i++) {
			coordinates.push_back(Point(coor[i*3], coor[i*3+1]));
			curvePos.push_back(coor[i*3+2]);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		//cout << "kończy procesor " << rank << endl;
	}
	vector <pdd> res;
	for(int i=0;i<pointCount;i++) {
		res.clear();
		query(Point(0.0, 0.0), 1.0, coordinates[i], radius, 1.0, res);
		odpowiedzi.push_back(pivpdd(i,res));	
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
	/*for(int i=0;i<(int)assignements.size();i++) {
		cout << "procesor " << rank << " ma sasiada od " << assignements[i].p << " do " << assignements[i].k << "ktory jest procesorem " << assignements[i].procesor << endl;
		assert(assignements[i].p-1e-12 <= assignements[i].k);
	}*/

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

	/*MPI_Finalize();
	return 0;*/

	


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
	
	MPI_Request *kupaw = new MPI_Request[size];

	for(int i=0;i<size;i++) if(i!=rank && wyn[i].size() > 0) {
		double *buf = new double[wyn[i].size()*3];
		for(int k=0;k<(int)wyn[i].size();k++) {
			buf[k*3] = wyn[i][k].x;
			buf[k*3+1] = wyn[i][k].y;
			buf[k*3+2] = wyn77[i][k];
		}
		MPI_Isend(buf, wyn[i].size()*3, MPI_DOUBLE, i,7,MPI_COMM_WORLD, &kupaw[i]);
	}
	
	MPI_Request *bagnov = new MPI_Request[size];

	vector<pdp> box;
	
	
	double **buff = new double*[size];

	for(int i=0;i<size;i++) if(i!=rank && ile_punktow_ma_dla_mnie_sasiad[i] > 0) {
		buff[i] = new double[ile_punktow_ma_dla_mnie_sasiad[i]*3];
		MPI_Irecv(buff[i], ile_punktow_ma_dla_mnie_sasiad[i]*3, MPI_DOUBLE,i,
		              7, MPI_COMM_WORLD, &bagnov[i]);
	}
	


	for(int i=0;i<size;i++) if(i!=rank && wyn[i].size() > 0)  {
		MPI_Wait(&kupaw[i],&stat[i]);
	}
	for(int i=0;i<size;i++) if(i!=rank && ile_punktow_ma_dla_mnie_sasiad[i] > 0) {
		MPI_Wait(&bagnov[i],&stat[i+rank]);
	}
	

	for(int i=0;i<size;i++) {
		for(int j=0;j<ile_punktow_ma_dla_mnie_sasiad[i]*3;j+=3) {
			box.push_back(pdp(buff[i][j+2],Point(buff[i][j],buff[i][j+1])));
		}
	}

	sort(box.begin(),box.end());

	/*for(int i=0;i<box.size();i++) {
		cout << "jestem " << rank << " Mam sasiada!  : " << box[i].second.x << " " << box[i].second.y << endl;
	}*/

	/* Van der Vaals Force computing
	vector <double> forcex(pointCount);
	vector <double> forcey(pointCount);

	for(int i=0;i<pointCount;i++) {
		for(int j=0;j<pointCount;j++) if(i!=j) {
			forcex[i] += (coordinates[j].x-coordinates[i].x)/
			coordinates[i].distance(coordinates[j]) * force_van_der_vaals(coordinates[i], coordinates[j]);
			forcey[i] += (coordinates[j].y-coordinates[i].y)/
			coordinates[i].distance(coordinates[j]) * force_van_der_vaals(coordinates[i], coordinates[j]);
		}
	}

	for(int i=0;i<pointCount;i++) {
		for(int j=0;j<(int)box.size();j++) if(coordinates[i].distance(box[j].second) <= r) {
			forcex[i] += (box[j].second.x-coordinates[i].x)/
			coordinates[i].distance(box[j].second) * force_van_der_vaals(coordinates[i], box[j].second);
			forcey[i] += (box[j].second.y-coordinates[i].y)/
			coordinates[i].distance(box[j].second) * force_van_der_vaals(coordinates[i], box[j].second);
		}
	}
	*/
	/*for(int i=0;i<pointCount;i++) {
		cout << "punkt o wspolrzednych (" << coordinates[i].x << ", " << coordinates[i].y << ") -> " << "f(x) = " << forcex[i] << ", f(y) = " << forcey[i] << endl;
	}*/

	
	cout << "Jestem procesem " << rank << " i kończę\n";
	MPI_Finalize();
	return 0;

}

