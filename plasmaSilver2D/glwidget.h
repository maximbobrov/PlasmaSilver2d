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
    void setField( double **arr, int N_X, int N_Y, double dx, double dy);


private:
     double **arr = nullptr;
     int N_X;
     int N_Y;
     double dx;
     double dy;


};

#endif // GLWIDGET_H
