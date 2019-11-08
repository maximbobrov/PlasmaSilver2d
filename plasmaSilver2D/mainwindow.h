#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "qcustomplot.h"
#include <QGridLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QLineEdit>
#include <QScrollBar>
#include <QLabel>
#include <QVector>
#include <QCheckBox>
#include "glwidget.h"

#include "simulationdata.h"
#include "simulationsolver.h"
#include "simulationtools.h"


class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();

    struct plotStruct{
        QVector<double> x, y;
        double* arr;
        QString name;
        int size;
        double scale;
        bool visible;
    };

    struct storeStruct{
        QVector<plotStruct> plots;
        double time;
    };

public slots:
    void initData();
    void updateData();
    void simulateData(bool);
    void drawDebug(bool);
    void replotGraph(int);
private:
    QWidget* m_widget;
    GLWidget *glWidget;
    QCustomPlot* m_customPlot;
    QGridLayout* m_grid;
    QVBoxLayout* m_vLayoutCheckBoxes;
    QHBoxLayout* m_hLayout;
    QLineEdit* m_textStartTime;
    QLineEdit* m_textEndTime;
    QLineEdit* m_textDeltaTime;
    QScrollBar* m_scrollBar;
    QScrollBar* m_scallingBar;
    QVector<QCheckBox*> m_checkBoxes;
    QPushButton* m_simulateButton;
    QPushButton* m_debugButton;
    simulationData* m_data;
    solverNe* m_sNe;
    solverEnergy* m_sEn;
    solverPhi* m_sPhi;
    QVector<solverHeavySpicies*> m_sHeavy;
    int m_numberHeavySpicies;
    simulationData::simulationField2d* m_fNe;
    simulationData::simulationField2d* m_fEnergy;
    simulationData::simulationField2d* m_fPhi;
    QVector<simulationData::simulationField2d*> m_fHeavy;

    QVector<plotStruct> m_plots;
    QVector<storeStruct> m_storage;
    double m_maxY, m_minY;
    bool m_animStopped;
    double m_startTime,m_endTime;
    double m_time;
private:
    void saveInStorage();
    void addPlot(double* arr,char* name, int size, double scale = 1.0);
    void addPlotXY(double *arr,double*xx, char *name, int size, double scale = 1.0);
};

#endif // MAINWINDOW_H
