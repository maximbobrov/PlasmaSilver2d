#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>
#include <QMouseEvent>
#include <QEvent>

typedef struct {
  float x,y,z;
} XYZ;

class GLWidget : public QGLWidget
{
    Q_OBJECT
public:
    explicit GLWidget(QWidget *parent = 0);
    

    void initializeGL();
    void paintGL();
    void resizeGL(int w, int h);
    XYZ get_color(double val,double min,double max);
    void setField( double **arr, double **eps, int N_X, int N_Y, double dx, double dy, double sc);


private:
     double **arr = nullptr;
     double **eps = nullptr;
     int N_X;
     int N_Y;
     double dx;
     double dy;
     double sc;
};

#endif // GLWIDGET_H
