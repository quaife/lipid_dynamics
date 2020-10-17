// iagm_ves.h

#include <Eigen/Core>
#include <Curve.h>
#include <InterferenceVolume.h>

using namespace std;
using namespace Eigen;



class Iagm_ves : public IAGM::Curve
{
public:
	size_t m_qsize;
	size_t m_psize;
	double *m_xStart;
	double *m_xEnd;
	double *m_xCtrl;
	size_t *m_xIndex;
        //double *m_xCtrl_copy;
	unsigned int *m_edges;        
	double *m_volumeGra;
	double *m_ids;
	double *m_vols;
	double m_volume;
	//size_t m_nv;
	//size_t m_nb;
	//size_t m_Npv;
	//size_t m_Npb;
        // constructor
        // n, number of points on each vesicle
        // m. number of vesicles
        // x1, vesicle start position
        // x2, vesicle end position
        // min_sep, minimal separation distance
	Iagm_ves(double n[], size_t m, double x1[], double x2[],double min_sep,unsigned int nbriters, size_t np, double tol, size_t nv, size_t nb, size_t Npv, size_t Npb,size_t nexten, double cellSize);
	~Iagm_ves();
	// update final positions
	void updateQ();
	void getVolume();
	void getVolumeGradient();
        // Derived Methods
	void updateSlaveMesh();
        void getSubspaceVector(double *in, double *out);
        size_t getNbrVertices() { return m_qsize; }
        size_t getNbrCtrlVertices() { return m_psize; }
        void getStartConfiguration(double *&x);
        void getEndConfiguration(double *&x);
        void getControlConfiguration(double *&x);
	void getIndex(size_t *&x);
        void copyControlConfiguration();
        void getMovableControlVertices(std::set<size_t> &movable);
        void getEdgeIndices(size_t &nbr, unsigned int *&e);
};
