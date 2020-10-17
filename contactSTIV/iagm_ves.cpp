// iagm_ves.cpp
//
#include <iostream>
#include "iagm_ves.h"
//using namespace std;
Iagm_ves::Iagm_ves(double n[], size_t m, double x1[], double x2[], double min_sep, unsigned int nbriters, size_t np, double tol, size_t nv, size_t nb, size_t Npv, size_t Npb, size_t nexten, double cellSize)
{

	m_qsize = np;
	m_psize = np;
        
	m_xStart = new double[2 * m_qsize];
	m_xEnd   = new double[2 * m_qsize];
	m_xCtrl  = new double[2 * m_psize];
	m_xIndex = new size_t[2 * m_qsize];


	m_volumeGra = new double[2 * m_qsize];
	m_ids = new double[2 * m_qsize];
	m_vols = new double[2 * m_qsize];

	m_edges = new unsigned int[2 * m_qsize];
	m_nv = nv;
	m_nb = nb;
	m_Npv = Npv;
	m_Npb = Npb;
	m_nexten = nexten;
	size_t nCount = 0;
	// Fill out edge data structure and copy x1, x2
	for (size_t i=0, currPt=0; i<m; i++)
	{
		size_t nn = (size_t)n[i];
		for(size_t j=0; j<nn-1; j++, currPt++)
		{
			m_edges[2*currPt  ] = currPt;
			m_edges[2*currPt+1] = currPt + 1;
                        
			m_xStart[2*currPt  ] = x1[2*nCount+j];
			m_xStart[2*currPt+1] = x1[2*nCount+j+nn];
                
			m_xCtrl[2*currPt  ] = x2[2*nCount+j];
			m_xCtrl[2*currPt+1] = x2[2*nCount+j+nn];

			m_xIndex[2*currPt  ] = i;
			m_xIndex[2*currPt+1] = i;
						
			m_volumeGra[2*currPt  ] = 0;
			m_volumeGra[2*currPt+1] = 0;
						
			m_ids[2*currPt  ] = 0;
			m_ids[2*currPt+1] = 0;
						
			m_vols[2*currPt  ] = 0;
			m_vols[2*currPt+1] = 0;

		}
                
		m_edges[2*currPt  ] = currPt;
		m_edges[2*currPt+1] = currPt+1-nn;

		m_xStart[2*currPt  ] = x1[2*nCount+nn-1];
		m_xStart[2*currPt+1] = x1[2*nCount+nn-1+nn];
                
		m_xCtrl[2*currPt  ] = x2[2*nCount+nn-1];
		m_xCtrl[2*currPt+1] = x2[2*nCount+nn-1+nn];
		
		m_xIndex[2*currPt  ] = i;
		m_xIndex[2*currPt+1] = i;
				
		m_volumeGra[2*currPt  ] = 0;
		m_volumeGra[2*currPt+1] = 0;
				
		m_ids[2*currPt  ] = 0;
		m_ids[2*currPt+1] = 0;
				
		m_vols[2*currPt  ] = 0;
		m_vols[2*currPt+1] = 0;
		
		currPt++;
		
		nCount = nCount + nn;
	}


        setCellSize(cellSize);
	setMinimumSeparation(min_sep);
	setMaxNbrIterations(nbriters);
	setTol(tol);

}

void Iagm_ves::updateQ()
{
//        std::cout<<getMinimumSeparation()<<std::endl;
        updateSlaveMesh();

        findAndRemoveInterference();

} 

void Iagm_ves::getVolume()
{
	updateSlaveMesh();
	IAGM::Intersections is;
	m_volume = 0;
	if(getCurveIntersections(is))
	{
		IAGM::InterferenceVolumes ivs;
		modelInterference(is, ivs);
		//std::cout<<"size of ivs "<<ivs.size()<<std::endl;
		set<size_t> movable;
		getMovableControlVertices(movable);
		int countId = 0;
		for(IAGM::InterferenceVolumesIterator ivItr=ivItr=ivs.begin(); ivItr!=ivs.end(); ++ivItr )
		{
			if(ivItr->getV()>=0)
				continue;
			countId++;
			ivItr->computeHandleGradient(this, movable);
			//std::cout<<"volume: "<<ivItr->getV()<<std::endl;    
			m_volume += ivItr->getV();
			for (set<size_t>::iterator sItr=movable.begin();sItr!=movable.end();++sItr)
			{
				for (size_t i=0; i<2; ++i)
				{
					m_volumeGra[2*(*sItr)+i] += ivItr->getHandleGradients()[2*(*sItr)+i];
					if(ivItr->getHandleGradients()[2*(*sItr)+0]!=0 ||  ivItr->getHandleGradients()[2*(*sItr)+1]!=0)
					{
						m_ids[2*(*sItr)+i] = countId;
						m_vols[2*(*sItr)+i] = ivItr->getV();
					}
//					std::cout<<ivItr->getHandleGradients()[2*(*sItr)+i]<<std::endl;
				}
			}
		
		}
		
	}
	std::cout<<"volume is "<<m_volume<<endl;
}

void Iagm_ves::updateSlaveMesh()
{

        for (size_t i=0;i<m_qsize;i++)
        {
                m_xEnd[2*i  ] = m_xCtrl[2*i  ];
                m_xEnd[2*i+1] = m_xCtrl[2*i+1];
        }
}


void Iagm_ves::getSubspaceVector(double *in, double *out)
{
        for (size_t i=0;i<m_qsize;i++)
        {
                out[2*i+0] = in[2*i+0];
                out[2*i+1] = in[2*i+1];
        }       
}


void Iagm_ves::getStartConfiguration(double *&x)
{
        x = m_xStart;
}



void Iagm_ves::getEndConfiguration(double *&x)
{
        x = m_xEnd;
}

void Iagm_ves::getControlConfiguration(double *&x)
{
        x = m_xCtrl;
}

void Iagm_ves::getIndex(size_t *&x)
{
	x = m_xIndex;
}

void Iagm_ves::copyControlConfiguration()
{
        return;
}

// need change when add boundary
void Iagm_ves::getMovableControlVertices(set<size_t> &movable)
{
        movable.clear();
        for (size_t i=0;i<m_psize;i++)
                movable.insert(i);
}

void Iagm_ves::getEdgeIndices(size_t &nbr, unsigned int *&e)
{
        e = m_edges;
        nbr = m_qsize;
}

Iagm_ves::~Iagm_ves()
{
	delete m_xStart;
	delete m_xEnd;
	delete m_xCtrl;
	delete m_xIndex;
	delete m_volumeGra;
	delete m_ids;
	delete m_vols;
	delete m_edges;
}
