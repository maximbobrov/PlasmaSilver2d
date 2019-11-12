#include "glwidget.h"
#include <math.h>
#include <QDebug>

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

XYZ GLWidget::get_color(double gval, double min, double max)
{
    const int nn=4;
    int i;
    double val;
    val=gval;
    if (val>max) val=max-0.00001;
    if (val<min) val=min+0.00001;

    XYZ col_table[5];
    XYZ res;

    col_table[0].x=0.0; col_table[0].y=0.0; col_table[0].z=1.0;
    col_table[1].x=0.0; col_table[1].y=1.0; col_table[1].z=1.0;
    col_table[2].x=0.0; col_table[2].y=1.0; col_table[2].z=0.0;
    col_table[3].x=1.0; col_table[3].y=1.0; col_table[3].z=0.0;
    col_table[4].x=1.0; col_table[4].y=0.0; col_table[4].z=0.0;

    double alpha;
    if ((max-min)>0.00001)
    {
        alpha=(val-min)/(max-min)*nn;
        i=(int)(alpha);
        alpha=alpha-i;
    }else
    {
        alpha=0.0;
        i=2;
    }
    res.x=col_table[i].x*(1-alpha)+col_table[i+1].x*alpha;
    res.y=col_table[i].y*(1-alpha)+col_table[i+1].y*alpha;
    res.z=col_table[i].z*(1-alpha)+col_table[i+1].z*alpha;


    glColor3f(res.x,res.y,res.z);
    return res;

}

void GLWidget::setField(double **iArr, double **iEps, int iN_X, int iN_Y, double idx, double idy, double isc = 1.0)
{
    arr = iArr;
    eps = iEps;
    N_X = iN_X;
    N_Y = iN_Y;
    dx =  idx;
    dy =  idy;
    sc = isc;
}

void GLWidget::paintGL()
{
    glColor3f(1,1,1);

    glClearColor (0.0, 0.0, 0.0, 0.0);
    glEnable(GL_LINE_SMOOTH);


    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();
    //glPushMatrix();
    glScalef(4,4,4);
    glDisable(GL_FOG);


    //glPopMatrix();
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_FOG);
    //glDisable(GL_LINE_SMOOTH);
    glDisable(GL_LIGHTING);

    //glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);

    /*glBegin(GL_TRIANGLE_STRIP);
    for (int i=0;i<=20;i++)
    {
        get_color(i/20.0,0.0,1.0);

        glVertex2f(-0.75+i*1.5/20,-0.94);
        glVertex2f(-0.75+i*1.5/20,-0.98);
    }
    glEnd();*/

    if(arr!=nullptr)
    {

        double l_2;
        double A = 1.0;
        glColor3f(1,1,1);
        for (int i=1;i<N_X-1;i++)
        {
            glBegin(GL_TRIANGLE_STRIP);
            for (int j=0;j<N_Y;j++)
            {
                l_2=sc*(arr[i][j]/A);
                glColor3f(l_2,l_2,-l_2);

                glVertex2f(dx*(i-N_X/2),dy*(j-N_Y/2));

                l_2=sc*(arr[i+1][j]/A);
                glColor3f(l_2,l_2,-l_2);
                glVertex2f(dx*(i+1-N_X/2),dy*(j-N_Y/2));
            }
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
        }

        glColor3f(1,1,1);

        glBegin(GL_LINE_LOOP);

        glVertex3f(dx*(-N_X/2),dy*(-N_Y/2),0);
        glVertex3f(dx*(N_X-1-N_X/2),dy*(-N_Y/2),0);
        glVertex3f(dx*(N_X-1-N_X/2),dy*(N_Y-1-N_Y/2),0);
        glVertex3f(dx*(-N_X/2),dy*(N_Y-1-N_Y/2),0);
        glEnd();
    }




    glEnable(GL_DEPTH_TEST);
    //glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);


}

void GLWidget::resizeGL(int w, int h)
{


    glClearColor (0.0, 0.0, 0.0, 0.0);

    /* set fill color to white */
    glColor3f(1.0, 1.0, 1.0);

    glViewport(0,0,w,h);

    /* set up standard orthogonal view with clipping */
    /* box as cube of side 2 centered at origin */
    /* This is default view and these statement could be removed */
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
