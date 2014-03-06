/*
  FLUIDS v.1 - SPH Fluid Simulator for CPU and GPU
  Copyright (C) 2008. Rama Hoetzlein, http://www.rchoetzlein.com

  ZLib license
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
*/

#include "gl_helper.h"

#include "point_set.h"

int PointSet::m_pcurr = -1;


PointSet::PointSet ()
{	
	m_GridRes.Set ( 0, 0, 0 );
	//myvertices;
	m_pcurr = -1;
	Reset ();
}

int PointSet::GetGridCell ( int x, int y, int z )
{
	return (int) ( (z*m_GridRes.y + y)*m_GridRes.x + x);
}

Point* PointSet::firstGridParticle ( int gc, int& p )
{
	m_pcurr = m_Grid [ gc ];
	if ( m_pcurr == -1 ) return 0x0;
	p = m_pcurr;
	return (Point*) (mBuf[0].data + m_pcurr * mBuf[0].stride);
}

Point* PointSet::nextGridParticle ( int& p )
{
	Point* pnt = 0x0;
	if ( m_pcurr != -1 ) {
		pnt = (Point*) (mBuf[0].data + m_pcurr * mBuf[0].stride);
		p = m_pcurr;
		m_pcurr = pnt->next;
	}	
	return pnt;
}

unsigned short* PointSet::getNeighborTable ( int n, int& cnt )
{
	cnt = m_NC[n];
	if ( cnt == 0 ) return 0x0;
	return &m_Neighbor[n][0];
}

float PointSet::GetValue ( float x, float y, float z )
{
	float dx, dy, dz, dsq;
	float sum;
	int pndx;
	Point* pcurr;
	float R2 = 1.8*1.8;

	Grid_FindCells ( Vector3DF(x,y,z), m_GridCellsize/2.0 );

	int cnt = 0;
	sum = 0.0;
	for (int cell=0; cell < 8; cell++ ) {
		if ( m_GridCell[cell] != -1 ) {
			pndx = m_Grid [ m_GridCell[cell] ];
			while ( pndx != -1 ) {					
				pcurr = (Point*) (mBuf[0].data + pndx*mBuf[0].stride);
				dx = x - pcurr->pos.x;
				dy = y - pcurr->pos.y;
				dz = z - pcurr->pos.z;
				dsq = dx*dx+dy*dy+dz*dz;		
				if ( dsq < R2 ) sum += R2 / dsq;
				pndx = pcurr->next;
			}	
		}
	}
	return sum;	
}	
Vector3DF PointSet::GetGradient ( float x, float y, float z )
{
	Vector3DF norm;  
	float dx, dy, dz, dsq;
	float sum;
	int pndx;
	Point* pcurr;
	float R2 = (m_GridCellsize/2.0)*(m_GridCellsize/2.0);

	Grid_FindCells ( Vector3DF(x,y,z), m_GridCellsize/2.0 );

	int cnt = 0;
	sum = 0.0;
	norm.Set (0,0,0);
	for (int cell=0; cell < 8; cell++ ) {
		if ( m_GridCell[cell] != -1 ) {
			pndx = m_Grid [ m_GridCell[cell] ];
			while ( pndx != -1 ) {					
				pcurr = (Point*) (mBuf[0].data + pndx*mBuf[0].stride);
				dx = x - pcurr->pos.x;
				dy = y - pcurr->pos.y;
				dz = z - pcurr->pos.z;
				dsq = dx*dx+dy*dy+dz*dz;				
				if ( dsq > 0 && dsq < R2 ) {
					dsq = 2.0*R2 / (dsq*dsq);
					norm.x += dx * dsq;
					norm.y += dy * dsq;
					norm.z += dz * dsq;						
				}
				pndx = pcurr->next;
			}	
		}
	}
	norm.Normalize ();	
	return norm;
}

DWORD PointSet::GetColor ( float x, float y, float z )
{
	Vector3DF clr;  
	float dx, dy, dz, dsq;
	float sum;
	int pndx;
	Point* pcurr;
	float R2 = (m_GridCellsize/2.0)*(m_GridCellsize/2.0);

	Grid_FindCells ( Vector3DF(x,y,z), m_GridCellsize/2.0 );

	int cnt = 0;
	sum = 0.0;
	clr.Set (0,0,0);
	for (int cell=0; cell < 8; cell++ ) {
		if ( m_GridCell[cell] != -1 ) {
			pndx = m_Grid [ m_GridCell[cell] ];
			while ( pndx != -1 ) {					
				pcurr = (Point*) (mBuf[0].data + pndx*mBuf[0].stride);
				dx = x - pcurr->pos.x;
				dy = y - pcurr->pos.y;
				dz = z - pcurr->pos.z;
				dsq = dx*dx+dy*dy+dz*dz;				
				if ( dsq < R2 ) {
					dsq = 2.0*R2 / (dsq*dsq);					
					clr.x += RED(pcurr->clr) * dsq;
					clr.y += GRN(pcurr->clr) * dsq;
					clr.z += BLUE(pcurr->clr) * dsq;						
				}
				pndx = pcurr->next;
			}	
		}
	}
	clr.Normalize ();
	return COLORA(clr.x, clr.y, clr.z, 1.0);
}

void PointSet::Reset ()
{
	// Reset number of particles
//	ResetBuffer ( 0 );	

	m_Time = 0;
	m_DT = 0.1;
	m_Param[POINT_GRAV] = 100.0;
	m_Param[PLANE_GRAV] = 0.0;
	
	m_Vec[ POINT_GRAV_POS].Set(0,0,50.0);
	m_Vec[ PLANE_GRAV_DIR].Set(0,0,-9.8);
	m_Vec[ EMIT_RATE ].Set ( 1, 10, 0 );
	m_Vec[ EMIT_POS ].Set ( 50, 0, 35 );	
	m_Vec[ EMIT_ANG ].Set ( 90, 45, 50.0 );	
	m_Vec[ EMIT_DANG ].Set ( 0, 0, 0 );	
	m_Vec[ EMIT_SPREAD ].Set ( 4, 4, 1 );	
}

void PointSet::Initialize ( int mode, int total )
{
	switch (mode) {
	case BPOINT: {
		FreeBuffers ();
		AddBuffer ( BPOINT, sizeof ( Point ), total );
		AddAttribute ( 0, "pos", sizeof ( Vector3DF ), false );
		AddAttribute ( 0, "color", sizeof ( DWORD ), false );
		AddAttribute ( 0, "type", sizeof ( unsigned short), false );
		Reset ();
		} break;
	
	case BPARTICLE: {
		FreeBuffers ();
		AddBuffer ( BPARTICLE, sizeof ( Particle ), total );
		AddAttribute ( 0, "pos", sizeof ( Vector3DF ), false );	
		AddAttribute ( 0, "color", sizeof ( DWORD ), false );
		AddAttribute ( 0, "vel", sizeof ( Vector3DF ), false );
		AddAttribute ( 0, "ndx", sizeof ( unsigned short ), false );
		AddAttribute ( 0, "age", sizeof ( unsigned short ), false );
		AddAttribute ( 0, "type", sizeof ( unsigned short), false );
		Reset ();
		} break;
	}


}

int PointSet::AddPoint ()
{
	xref ndx;	
	AddElem ( 0, ndx );	
	return ndx;
}

int PointSet::AddPointReuse ()
{
	xref ndx;	
	if ( NumPoints() < mBuf[0].max-1 )
		AddElem ( 0, ndx );
	else
		RandomElem ( 0, ndx );
	return ndx;
}



int PointSet::AddVolume ( Vector3DF min, Vector3DF max, float spacing )
{
	Vector3DF pos;
	Point* p;	
	float dx, dy, dz;
	dx = max.x-min.x;
	dy = max.y-min.y;
	dz = max.z-min.z;
	int cnt = 0;
	for (float z = max.z; z >= min.z; z -= spacing ) {
		for (float y = min.y; y <= max.y; y += spacing ) {	
			for (float x = min.x; x <= max.x; x += spacing ) {
				p = GetPoint ( AddPointReuse () );
				pos.Set ( x, y, z);
				//pos.x += -0.05 + float( rand() * 0.1 ) / RAND_MAX;
				//pos.y += -0.05 + float( rand() * 0.1 ) / RAND_MAX;
				//pos.z += -0.05 + float( rand() * 0.1 ) / RAND_MAX;
				p->pos = pos;	
				p->type = 0;
				p->clr = COLORA( (x-min.x)/dx, (y-min.y)/dy, (z-min.z)/dz, 1);
				cnt++;
			}
		}
	}
	sphPoints = cnt ;
	return sphPoints;
}

void PointSet::AddControlVolume ( std::vector<glm::vec3> *allvertices )
{
    Vector3DF pos;
	Point* p;	
	int dimx=5;
	int dimy=5;
	float dimz=10;
	int ccnt=0;
	//vec3 a = myvertices. ; 
	glm::vec3 startpoint(0,0.0,-5.0);
	glm::vec3 direction(0,0,1.0);
	glm::vec3 vertex1(0,1.0,2.0);
	glm::vec3 vertex2(-1.0,-1.0,0);
	glm::vec3 vertex3(1.0,-1.0,0);
	glm::vec3 v1,v2,v3;
	int a = allvertices->size();
	double value=0;
	value=PointSet::Test_RayPolyIntersect(startpoint,direction,vertex1,vertex2,vertex3);
	int count=0;
	float var=1.0;
	//glm::vec3 *a = myvertices;
	for (float x = -5; x < 18; x += 1.0)
	{
		for (float y = -5; y <= 14; y += 1.0)
		{	
			 for (float z = -5; z <= 25; z += var )
			 {
				 /*if(z>=-5.0 && z <2.0)
					 var =1.0;
				  else if(z>=2.0 && z < 8.0)
					 var =0.70;
				   else if(z>=8.0 && z < 15.0)
					 var =0.6;
				    else if(z >= 15.0 && z <=25.0)
					 var =0.3;*/

				 count = 0 ;
				 startpoint=glm::vec3((float)x,(float)y,(float)z);

				 direction=glm::vec3(1.0,0.1,0);
				 for(int v=0 ;v < a ; v+=3)
				 {
					 //checking each vertex in the grid with all the 192 triangles of mesh
					 v1=(*allvertices)[v];
					 v2=(*allvertices)[v+1];
					 v3=(*allvertices)[v+2];
					 vertex1= vec3(v1[0],v1[1],v1[2]);
					 vertex2= vec3(v2[0],v2[1],v2[2]);
					 vertex3= vec3(v3[0],v3[1],v3[2]);
					 value=PointSet::Test_RayPolyIntersect(startpoint,direction,vertex1,vertex2,vertex3);
						if(value != -1.0)
						{
							 count++;
						}

				 }

				 //odd
				 if( count % 2 == 1)
				 {
							 //inside the mesh
						p = GetPoint ( AddPointReuse () );
						pos.Set ( x,y,z+10.0f);
						p->pos = pos;	
						p->type = 1;
						p->clr = COLORA( 1.0,0.3,0.3, 0);  ///   0.1,0.3,1.0,1.0 
						
						ccnt++;
				 }
				 //var=var-0.1;
			}
			// var = 1.0;
			 
		}

	}	

}

void PointSet::Draw ( float* view_mat, float rad )
{
	char* dat;
	Point* p;
	glEnable ( GL_NORMALIZE );	

	if ( m_Param[PNT_DRAWMODE] == 0 ) {
		glLoadMatrixf ( view_mat );
		dat = mBuf[0].data;	
		for (int n = 0; n < NumPoints(); n++) {
			p = (Point*) dat;
			glPushMatrix ();
			glTranslatef ( p->pos.x, p->pos.y, p->pos.z );		
			glScalef ( 0.2, 0.2, 0.2 );	
			if(p->type == 0)
			glColor4f ( 0.1,0.3,1.0,1.0 );//glColor4f ( RED(p->clr), GRN(p->clr), BLUE(p->clr), ALPH(p->clr) );
			else
            glColor4f ( RED(p->clr), GRN(p->clr), BLUE(p->clr), 0.0 );
			drawSphere ();
			glPopMatrix ();		
			dat += mBuf[0].stride;
		}	
	} else if ( m_Param[PNT_DRAWMODE] == 1 ) {
		glLoadMatrixf ( view_mat );
		dat = mBuf[0].data;
		glBegin ( GL_POINTS );
		for (int n=0; n < NumPoints(); n++) {
			p = (Point*) dat;
			glColor3f ( RED(p->clr), GRN(p->clr), BLUE(p->clr) );			
			glVertex3f ( p->pos.x, p->pos.y, p->pos.z );			
			dat += mBuf[0].stride;
		}
		glEnd ();
	}
}

void PointSet::Emit ( float spacing )
{
	Particle* p;
	Vector3DF dir;
	Vector3DF pos;
	float ang_rand, tilt_rand;
	float rnd = m_Vec[EMIT_RATE].y * 0.15;	
	int x = (int) sqrt(m_Vec[EMIT_RATE].y);

	for ( int n = 0; n < m_Vec[EMIT_RATE].y; n++ ) {
		ang_rand = (float(rand()*2.0/RAND_MAX) - 1.0) * m_Vec[EMIT_SPREAD].x;
		tilt_rand = (float(rand()*2.0/RAND_MAX) - 1.0) * m_Vec[EMIT_SPREAD].y;
		dir.x = cos ( ( m_Vec[EMIT_ANG].x + ang_rand) * DEGtoRAD ) * sin( ( m_Vec[EMIT_ANG].y + tilt_rand) * DEGtoRAD ) * m_Vec[EMIT_ANG].z;
		dir.y = sin ( ( m_Vec[EMIT_ANG].x + ang_rand) * DEGtoRAD ) * sin( ( m_Vec[EMIT_ANG].y + tilt_rand) * DEGtoRAD ) * m_Vec[EMIT_ANG].z;
		dir.z = cos ( ( m_Vec[EMIT_ANG].y + tilt_rand) * DEGtoRAD ) * m_Vec[EMIT_ANG].z;
		pos = m_Vec[EMIT_POS];
		pos.x += spacing * (n/x);
		pos.y += spacing * (n%x);
		
		p = (Particle*) GetElem( 0, AddPointReuse () );
		p->pos = pos;
		p->vel = dir;
		p->vel_eval = dir;
		p->age = 0;	
		p->type = 0;
		p->clr = COLORA ( m_Time/10.0, m_Time/5.0, m_Time /4.0, 1 );
	}
}


void PointSet::Run ()
{
	if ( m_Vec[EMIT_RATE].x > 0 && ++m_Frame >= (int) m_Vec[EMIT_RATE].x ) {
		m_Frame = 0;
		Emit ( 1.0 ); 
	}
	Advance();
}

void PointSet::Advance ()
{
	char* dat;
	Particle* p;
	Vector3DF vnext, accel, norm;

	vnext = m_Vec[EMIT_DANG];
	vnext *= m_DT;
	m_Vec[EMIT_ANG] += vnext;
		
	dat = mBuf[0].data;		
	for ( int c = 0; c < NumPoints(); c++ ) {		
		p = (Particle*) dat;

		accel.Set (0, 0, 0);

		// Plane gravity
		if ( m_Param[PLANE_GRAV] > 0) 
			accel += m_Vec[PLANE_GRAV_DIR];

		// Point gravity
		if ( m_Param[POINT_GRAV] > 0 ) {
			norm.x = ( p->pos.x - m_Vec[POINT_GRAV_POS].x );
			norm.y = ( p->pos.y - m_Vec[POINT_GRAV_POS].y );
			norm.z = ( p->pos.z - m_Vec[POINT_GRAV_POS].z );
			norm.Normalize ();
			norm *= m_Param[POINT_GRAV];
			accel -= norm;
		}
		if(p->type == 0){
		// Leapfrog Integration ----------------------------
		vnext = accel;		
		vnext *= m_DT;
		vnext += p->vel;								// v(t+1/2) = v(t-1/2) + a(t) dt
		p->vel_eval = p->vel;
		p->vel_eval += vnext;
		p->vel_eval *= 0.5;								// v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5		used to compute forces later
		p->vel = vnext;
		vnext *= m_DT;
		p->pos += vnext;								// p(t+1) = p(t) + v(t+1/2) dt
		}
		// Euler integration -------------------------------
		// accel += m_Gravity;
		// accel *= m_DT;
		// mParticles[c].vel += accel;				// v(t+1) = v(t) + a(t) dt
		// mParticles[c].vel_eval += accel;
		// mParticles[c].vel_eval *= m_DT/d;
		// mParticles[c].pos += mParticles[c].vel_eval;
		// mParticles[c].vel_eval = mParticles[c].vel; 
		
		dat += mBuf[0].stride;
	}

	m_Time += m_DT;
}

// Ideal grid cell size (gs) = 2 * smoothing radius = 0.02*2 = 0.04
// Ideal domain size = k*gs/d = k*0.02*2/0.005 = k*8 = {8, 16, 24, 32, 40, 48, ..}
//    (k = number of cells, gs = cell size, d = simulation scale)
void PointSet::Grid_Setup ( Vector3DF min, Vector3DF max, float sim_scale, float cell_size, float border )
{
	float world_cellsize = cell_size / sim_scale;
	m_Grid.clear ();
	m_GridMin = min;	m_GridMin -= border;
	m_GridMax = max;	m_GridMax += border;
	m_GridSize = m_GridMax;
	m_GridSize -= m_GridMin;
	m_GridCellsize = world_cellsize;
	m_GridRes.x = ceil ( m_GridSize.x / world_cellsize );		// Determine grid resolution
	m_GridRes.y = ceil ( m_GridSize.y / world_cellsize );
	m_GridRes.z = ceil ( m_GridSize.z / world_cellsize );
	m_GridSize.x = m_GridRes.x * cell_size / sim_scale;				// Adjust grid size to multiple of cell size
	m_GridSize.y = m_GridRes.y * cell_size / sim_scale;
	m_GridSize.z = m_GridRes.z * cell_size / sim_scale;
	m_GridDelta = m_GridRes;		// delta = translate from world space to cell #
	m_GridDelta /= m_GridSize;
	m_GridTotal = (int)(m_GridSize.x * m_GridSize.y * m_GridSize.z);

	m_Grid.clear ();
	m_GridCnt.clear ();
	m_GridControlCnt.clear();

	m_Grid.reserve ( m_GridTotal );
	m_GridCnt.reserve ( m_GridTotal );
	m_GridControlCnt.reserve ( m_GridTotal );
	for (int n=0; n < m_GridTotal; n++) {
		m_Grid.push_back ( -1 );
		m_GridCnt.push_back ( 0 );
		m_GridControlCnt.push_back ( 0 );
	}

}

void PointSet::Grid_Draw ( float* view_mat )
{
	float clr;
	int cx, cy, cz;
	float x1, y1, z1;
	float x2, y2, z2;
	int g = 0;

	glLoadMatrixf ( view_mat );
	glColor3f ( 0.7, 0.7, 0.7 );

	glBegin ( GL_LINES );	




	cz = 0;
	for ( cz = 0; cz < m_GridRes.z; cz++ ) {
	for ( cy = 0; cy < m_GridRes.y; cy++ ) {
	for ( cx = 0; cx < m_GridRes.x; cx++ ) {

		// Cell is not empty. Process it.
		//if ( m_Grid[g] != 0x0 ) {
		//	clr = m_GridCnt[g]/30.0;
			clr = 0.25;
			if ( clr <0.25) clr =0.25;
			if ( clr >1) clr =1 ;
			glColor3f ( clr, clr, clr );
			x1 = (cx * m_GridDelta.x) + m_GridMin.x;		x2 = ((cx+1) * m_GridDelta.x) + m_GridMin.x;
			y1 = (cy * m_GridDelta.y) + m_GridMin.y;		y2 = ((cy+1) * m_GridDelta.y) + m_GridMin.y;
			z1 = (cz * m_GridDelta.z) + m_GridMin.z;		z2 = ((cz+1) * m_GridDelta.z) + m_GridMin.z;
			glVertex3f ( x1, y1, z1 );			glVertex3f ( x2, y1, z1 );
			glVertex3f ( x2, y1, z1 );			glVertex3f ( x2, y2, z1 );
			glVertex3f ( x2, y2, z1 );			glVertex3f ( x1, y2, z1 );
			glVertex3f ( x1, y2, z1 );			glVertex3f ( x1, y1, z1 );
			glVertex3f ( x1, y1, z2 );			glVertex3f ( x2, y1, z2 );
			glVertex3f ( x2, y1, z2 );			glVertex3f ( x2, y2, z2 );
			glVertex3f ( x2, y2, z2 );			glVertex3f ( x1, y2, z2 );
			glVertex3f ( x1, y2, z2 );			glVertex3f ( x1, y1, z2 );
			glVertex3f ( x1, y1, z1 );			glVertex3f ( x1, y1, z2 );
			glVertex3f ( x1, y2, z1 );			glVertex3f ( x1, y2, z2 );
			glVertex3f ( x2, y2, z1 );			glVertex3f ( x2, y2, z2 );
			glVertex3f ( x2, y1, z1 );			glVertex3f ( x2, y1, z2 );
		//}
		//g++;
	}
	}
	}

	glEnd ();
}

void PointSet::Grid_InsertParticles ()
{
	char *dat1, *dat1_end;
	Point *p;
	int gs;
	int gx, gy, gz;
	
	dat1_end = mBuf[0].data + NumPoints()*mBuf[0].stride;
	for ( dat1 = mBuf[0].data; dat1 < dat1_end; dat1 += mBuf[0].stride ) 
		((Point*) dat1)->next = -1;	
//initialize
	for (int n=0; n < m_GridTotal; n++)
	{
		m_Grid[n] = -1;
		m_GridCnt[n] = 0;
		m_GridControlCnt[n] = 0;
	}

	dat1_end = mBuf[0].data + NumPoints()*mBuf[0].stride;
	int n = 0;
	for ( dat1 = mBuf[0].data; dat1 < dat1_end; dat1 += mBuf[0].stride ) {
		p = (Point*) dat1;
		gx = (int)( (p->pos.x - m_GridMin.x) * m_GridDelta.x);		// Determine grid cell
		gy = (int)( (p->pos.y - m_GridMin.y) * m_GridDelta.y);
		gz = (int)( (p->pos.z - m_GridMin.z) * m_GridDelta.z);
		//grid index gs
		gs = (int)( (gz*m_GridRes.y + gy)*m_GridRes.x + gx);		
		if ( gs >= 0 && gs < m_GridTotal )
		{
			p->next = m_Grid[gs];// pointer to the previous particle in the list of the recent particle
			m_Grid[gs] = n; // stores most recent particle index of that particular cell
			//that particular particle maintains a list of all the previous particles in the cell
			//if(p->type == 0)
			//{
			////m_GridCnt[gs]++; //stores no of particles in one grid cell
			//}
			//
			//if(p->type == 1)
			//{
			//m_GridControlCnt[gs]++; //stores no of particles in one grid cell
			//}

		}
		n++;
	}
}

int PointSet::Grid_FindCell ( Vector3DF p )
{
	int gc;
	Vector3DI cell;
	cell.x = (int) (p.x - m_GridMin.x) * m_GridDelta.x;
	cell.y = (int) (p.y - m_GridMin.y) * m_GridDelta.y;
	cell.z = (int) (p.z - m_GridMin.z) * m_GridDelta.z;
	gc = (int)( (cell.z*m_GridRes.y + cell.y)*m_GridRes.x + cell.x);
	if ( gc < 0 || gc > m_GridTotal ) return -1;
	return gc;
}

void PointSet::Grid_FindCells ( Vector3DF p, float radius )
{
	Vector3DI sph_min;

	// Compute sphere range
	sph_min.x = (int)((-radius + p.x - m_GridMin.x) * m_GridDelta.x);
	sph_min.y = (int)((-radius + p.y - m_GridMin.y) * m_GridDelta.y);
	sph_min.z = (int)((-radius + p.z - m_GridMin.z) * m_GridDelta.z);
	if ( sph_min.x < 0 ) sph_min.x = 0;
	if ( sph_min.y < 0 ) sph_min.y = 0;
	if ( sph_min.z < 0 ) sph_min.z = 0;

	m_GridCell[0] = (int)((sph_min.z * m_GridRes.y + sph_min.y) * m_GridRes.x + sph_min.x);
	m_GridCell[1] = m_GridCell[0] + 1;
	m_GridCell[2] = (int)(m_GridCell[0] + m_GridRes.x);
	m_GridCell[3] = m_GridCell[2] + 1;

	if ( sph_min.z+1 < m_GridRes.z ) {
		m_GridCell[4] = (int)(m_GridCell[0] + m_GridRes.y*m_GridRes.x);
		m_GridCell[5] = m_GridCell[4] + 1;
		m_GridCell[6] = (int)(m_GridCell[4] + m_GridRes.x);
		m_GridCell[7] = m_GridCell[6] + 1;
	}
	if ( sph_min.x+1 >= m_GridRes.x ) {
		m_GridCell[1] = -1;		m_GridCell[3] = -1;		
		m_GridCell[5] = -1;		m_GridCell[7] = -1;
	}
	if ( sph_min.y+1 >= m_GridRes.y ) {
		m_GridCell[2] = -1;		m_GridCell[3] = -1;
		m_GridCell[6] = -1;		m_GridCell[7] = -1;
	}
}



// Does the ray triangle intersection test based on the 560 code 
double PointSet::Test_RayPolyIntersect(glm::vec3 const& P0, glm::vec3 const& V0, glm::vec3 const& p11,   glm::vec3 const& p21,glm::vec3 const& p31)
{


	
	vec3 P = vec3(P0[0],P0[1],P0[2]);

	
	vec3 V = vec3(V0[0],V0[1],V0[2]);
	
	//glm mat is written in coloumns 
	vec4 pa = vec4(p11,1);
	vec3 p1 = vec3(pa[0],pa[1],pa[2]);
	vec4 pb =  vec4(p21,1);
	vec3 p2 = vec3(pb[0],pb[1],pb[2]);
	vec4 pc =  vec4(p31,1);
	vec3 p3 = vec3(pc[0],pc[1],pc[2]);
	

	vec3 n ;
	n = normalize(cross((p2-p1),(p3-p1)));
	double thit ; 
	thit  =  (dot(p1,n) - dot(P,n))/(dot(V,n)) ;
	if(thit < 0 )
		 return -1 ;

	// check if the intersection with plane was inside triangle ,if not then output -1 
	float sa,ta ;
	vec3 r,w;
	r = P + (float) thit * V ;
	w = r - p1 ;
	// Now using the parametric representation of a plane and finding the s and t values of the equation 
	// V(s,t) = V0 + s * ( V1 - V0) + t * (V2 - V0);
	// Find the s and t values using the ray equation . If s >=0 , t >=0 & s+t <= 0 ,then the point lies inside the triangle
	vec3 u,v;
	u = p2-p1;
	v = p3-p1;
	float den = pow(dot(u,v),2) - (dot(u,u) * dot(v,v)) ; 
	sa = (dot(u,v)*dot(w,v)  -  dot(v,v) * dot(w,u))/den ;
	ta = (dot(u,v)*dot(w,u)  -  dot(u,u) * dot(w,v))/den ;

	if ((sa >= 0) && (ta >= 0) && (sa+ta <= 1) )
	return thit;
	else
	return -1;
}
