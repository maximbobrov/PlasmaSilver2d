#include "glwidget.h"
#include <math.h>
#include <QDebug>
#include <QString>

GLWidget::GLWidget(QWidget *parent) :
    QGLWidget(parent)
{

}
void GLWidget::initializeGL()
{
    glClearColor (0.0, 0.0, 0.0, 0.0);
    /* set up standard orthogonal view with clipping */
    /* box as cube of side 2 centered at origin */
    /* This is default view and these statement could be removed */
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    glOrtho(-1.0, 1.0, -1.0, 1.0, -15, 15);
    glMatrixMode (GL_MODELVIEW);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_FOG);
    glEnable(GL_COLOR_MATERIAL);

}

void GLWidget::draw_text(double x, double y, double z, QString txt)
{
    glDisable(GL_LIGHTING);
    glDisable(GL_FOG);
    glColor3f(1,1,1);

    renderText(x, y, z, txt, QFont("Arial", 10, 0, false) );
    glEnable(GL_LIGHTING);
}


XYZ GLWidget::get_color(double gval, double min, double max)
{
    const int nn=4;
    int i;
    double val;
    val=gval;
    if (val>max) val=max;
    if (val<min) val=min;

    XYZ col_table[5];
    XYZ res;

    col_table[0].x = 0.0; col_table[0].y = 0.0; col_table[0].z = 1.0;
    col_table[1].x = 0.0; col_table[1].y = 1.0; col_table[1].z = 1.0;
    col_table[2].x = 0.0; col_table[2].y = 1.0; col_table[2].z = 0.0;
    col_table[3].x = 1.0; col_table[3].y = 1.0; col_table[3].z = 0.0;
    col_table[4].x = 1.0; col_table[4].y = 0.0; col_table[4].z = 0.0;

    double alpha;
    if ((max-min) > 1e-35)
    {
        alpha=(val-min)/(max-min)*nn;
        i=(int)(alpha);
        alpha=alpha-i;
    }
    else
    {
        alpha=0.0;
        i=2;
    }
    res.x = col_table[i].x * (1 - alpha) + col_table[i+1].x * alpha;
    res.y = col_table[i].y * (1 - alpha) + col_table[i+1].y * alpha;
    res.z = col_table[i].z * (1 - alpha) + col_table[i+1].z * alpha;

    glColor3f(res.x,res.y,res.z);
    return res;

}

void GLWidget::setField(double **iArr, double **iEps, double **iElectrod, int iN_X, int iN_Y, double idx, double idy, double iCrossNum, double isc = 1.0)
{
    arr = iArr;
    eps = iEps;
    electrod = iElectrod;
    N_X = iN_X;
    N_Y = iN_Y;
    dx =  idx;
    dy =  idy;
    sc = isc;
    crossNum = iCrossNum;
}
#include "mon/mcm.h"
void GLWidget::paintGL()
{
    glColor3f(1,1,1);

    glClearColor (0.0, 0.0, 0.0, 0.0);
    glEnable(GL_LINE_SMOOTH);


    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();
    glScalef(4,4,4);
    glDisable(GL_FOG);


    glDisable(GL_DEPTH_TEST);
    glDisable(GL_FOG);
    glDisable(GL_LIGHTING);

    int n=20;

    /*if(arr!=nullptr)
    {
        double max = 1e-35;
        double min = 1e35;
        glColor3f(1,1,1);

        for (int i=0;i<N_X;i++)
        {
            for (int j=0;j<N_Y;j++)
            {
                max = arr[i][j] > max ? arr[i][j] : max;
                min = arr[i][j] < min ? arr[i][j] : min;
            }
        }

        glBegin(GL_LINES);
        glColor3f(1,0,0);
        glLineWidth(3);
        glVertex2f(dx*(-N_X/2-3),dy*(crossNum-N_Y/2));
        glVertex2f(dx*(-N_X/2+N_X-1+3),dy*(crossNum-N_Y/2));
        glEnd();


        for (int i=0;i<N_X-1;i++)
        {
            glBegin(GL_TRIANGLE_STRIP);
            for (int j=0;j<N_Y;j++)
            {
                get_color(sc*arr[i][j],min,max);
                glVertex2f(dx*(i-N_X/2),dy*(j-N_Y/2));

                get_color(sc*arr[i+1][j],min,max);
                glVertex2f(dx*(i+1-N_X/2),dy*(j-N_Y/2));
            }
            //min = 0.0;
            //max = 1e11;
            glEnd();
            glPointSize(1);
            glBegin(GL_POINTS);
            glColor3f(1,1,1);

            for (int j=0;j<N_Y;j++)
            {
                if(eps[i][j] > 1.0)
                {
                    glVertex2f(dx*(i-N_X/2),dy*(j-N_Y/2));
                }
            }
            glEnd();

            glPointSize(1);
            glBegin(GL_POINTS);
            glColor3f(0.0,0.0,0.0);
            for (int j=0;j<N_Y;j++)
            {
                if(electrod[i][j] < 1.0)
                {
                    glVertex2f(dx*(i-N_X/2),dy*(j-N_Y/2));
                }
            }
            glEnd();



            glBegin(GL_TRIANGLE_STRIP);

            for (int i=0;i<=n;i++)
            {
                double l_3 = sc*((min + i*(max - min))/(n-1))/(max - min+ 1e-20 );
                get_color(sc*(min + i*(max - min))/(n-1),min,max);
                glVertex2f((i-n/2)*0.5/20,0.0-0.2*0.94);
                glVertex2f((i-n/2)*0.5/20,0.0-0.2*0.98);
            }
            glEnd();



        }

        glColor3f(1,1,1);

        glBegin(GL_LINE_LOOP);

        glVertex3f(dx*(-N_X/2),dy*(-N_Y/2),0);
        glVertex3f(dx*(N_X-1-N_X/2),dy*(-N_Y/2),0);
        glVertex3f(dx*(N_X-1-N_X/2),dy*(N_Y-1-N_Y/2),0);
        glVertex3f(dx*(-N_X/2),dy*(N_Y-1-N_Y/2),0);
        glEnd();

        draw_text((-n/2)*0.5/20,0.0-0.2*0.9, 0.0,  QString().number(min));
        draw_text((+n/2)*0.5/20,0.0-0.2*0.9, 0.0,  QString().number(max));
    }*/

    /**********Monte**********/
    using namespace monte;
    glLineWidth(1.0);
    glPointSize(2.0);


    /*glColor3f(0.1,0.1,0.1);
    glBegin(GL_QUADS);
    glVertex3f(0,0,0);
    glVertex3f(monte::dx*ncx,0,0);
    glVertex3f(monte::dx*ncx,dy*ncy,0);
    glVertex3f(0,dy*ncy,0);
    glEnd();*/


    printf("nsp=%d \n",nsp);
    printf("np[0]=%d \n",np[0]);
    printf("np[1]=%d \n",np[1]);
    glColor3f(1,0,0);


    double scc = 0.0000001;//0.0002;


    TWOD_One_2_One(sp_n_sm,   ncx, ncy);
    TWOD_One_2_One(sp_n_sm,   ncx, ncy);
    TWOD_One_2_One(sp_n_sm,   ncx, ncy);
    TWOD_One_2_One(sp_n_sm,   ncx, ncy);
    TWOD_One_2_One(sp_n_sm,   ncx, ncy);
    TWOD_One_2_One(sp_n_sm,   ncx, ncy);
    TWOD_One_2_One(sp_n_sm,   ncx, ncy);
    TWOD_One_2_One(sp_n_sm,   ncx, ncy);



    double sizeScale = 10;

   /* for (int i=0;i<ncx-1;i+=1)
    {
        glBegin(GL_TRIANGLE_STRIP);
        for (int j=0;j<ncy;j+=1)
        {
            glColor3f(sc * scc*sp_n_sm[i][j],-sc *scc*sp_n_sm[i][j],0);
            glVertex3f(sizeScale*monte::dx * i  - dx * N_X/2 , sizeScale*monte::dy * j  -  dy *N_Y/2, 0);
            glColor3f(sc * scc*sp_n_sm[i+1][j],-sc * scc*sp_n_sm[i+1][j],0);
            glVertex3f(sizeScale*monte::dx * (i+1)  - dx * N_X/2 , sizeScale*monte::dy * j  -  dy *N_Y/2, 0);
        }
        glEnd();
    }*/


   scc = 0.1000001;//0.0002;
    for (int i=0;i<ncx-1;i+=1)
    {
        glBegin(GL_TRIANGLE_STRIP);
        for (int j=0;j<ncy;j+=1)
        {
            glColor3f(sc * scc*Py_[i][j],-sc *scc*Py_[i][j],0);
            glVertex3f(sizeScale*monte::dx * i  - dx * N_X/2 , sizeScale*monte::dy * j  -  dy *N_Y/2, 0);
            glColor3f(sc * scc*Py_[i+1][j],-sc * scc*Py_[i+1][j],0);
            glVertex3f(sizeScale*monte::dx * (i+1)  - dx * N_X/2 , sizeScale*monte::dy * j  -  dy *N_Y/2, 0);
        }
        glEnd();
    }

    /*for (int i=0;i<ncx-1;i+=1)
    {
        glBegin(GL_TRIANGLE_STRIP);
        for (int j=0;j<ncy;j+=1)
        {
            glColor3f(sc * scc*ey[i][j],-sc *scc*ey[i][j],0);
            glVertex3f(sizeScale*monte::dx * i  - dx * N_X/2 , sizeScale*monte::dy * j  -  dy *N_Y/2, 0);
            glColor3f(sc * scc*ey[i+1][j],-sc * scc*ey[i+1][j],0);
            glVertex3f(sizeScale*monte::dx * (i+1)  - dx * N_X/2 , sizeScale*monte::dy * j  -  dy *N_Y/2, 0);
        }
        glEnd();
    }
*/

   /* glBegin(GL_TRIANGLE_STRIP);
    for (int i=0;i<pz_solver->m_p_num;i++)
    {
        glColor3f(sc * pz_solver->m_p[i].p/0.26,-sc * pz_solver->m_p[i].p/0.26,0);
        glVertex2f(sizeScale *pz_solver->m_p[i].r.x - dx * N_X/2 ,sizeScale *pz_solver->m_p[i].r.y -sizeScale *pz_solver->m_p[i].dl*0.5-  dy *N_Y/2);
        glVertex2f(sizeScale *pz_solver->m_p[i].r.x - dx * N_X/2 ,sizeScale *pz_solver->m_p[i].r.y +sizeScale *pz_solver->m_p[i].dl*0.5-  dy *N_Y/2);
    }
    glEnd();*/



    glBegin(GL_POINTS);
     for (int i=0;i<np[1];i+=1)
       {
        // printf("x=%e y=%e vx=%e vy=%e \n",x[0][i]*dx,y[0][i],vx[0][i],vy[0][i]);
double v2=vx[1][i]*vx[1][i]+vy[1][i]*vy[1][i];
//printf("v2=%e \n",v2);
           glColor3f(0,0,1.0);
         glVertex3f(sizeScale*monte::x[1][i]*monte::dx  - dx * N_X/2 ,
                       sizeScale*monte::y[1][i]*monte::dy  -  dy *N_Y/2,
                      0 );
       }
    glEnd();

    glBegin(GL_POINTS);
    for (int i=0;i<np[0];i+=1)
    {
        // printf("x=%e y=%e vx=%e vy=%e \n",x[0][i]*dx,y[0][i],vx[0][i],vy[0][i]);
        double v2=vx[0][i]*vx[0][i]+vy[0][i]*vy[0][i];
        //printf("v2=%e \n",v2);
        glColor3f(0,1.0,0);
        glVertex3f(sizeScale*monte::x[0][i]*monte::dx - dx * N_X/2,
                sizeScale*monte::y[0][i]*monte::dy -  dy *N_Y/2,
                0 );
    }

    /**********Monte**********/

    glEnd();

    glDisable(GL_DEPTH_TEST);


}

void GLWidget::resizeGL(int w, int h)
{
    glClearColor (0.0, 0.0, 0.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);

    glViewport(0,0,w,h);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    glOrtho(-w*1.0/h, w*1.0/h, -1.0, 1.0, -150.5, 150.5);
    glMatrixMode (GL_MODELVIEW);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    glEnable(GL_LINE_SMOOTH);

    glEnable(GL_FOG);
    glEnable(GL_DEPTH_TEST);
}
