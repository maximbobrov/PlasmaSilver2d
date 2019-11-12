
#include "mainwindow.h"
#include "qcustomplot.h"
#include <QGridLayout>
#include <QPushButton>
#include <QtCore>
#include <QVector>
#include <QDebug>
#include <stdio.h>

#define NZ 400

#define NX 129
#define NY 129



MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    m_widget=new QWidget();
    m_widget->setFixedHeight(700);
    m_widget->setFixedWidth(900);

    m_widget->setObjectName("central");
    glWidget = new GLWidget;



    m_customPlot = new QCustomPlot(this);
    m_grid = new QGridLayout(this);
    m_vLayoutCheckBoxes = new QVBoxLayout(this);
    m_hLayout = new QHBoxLayout(this);
    m_textStartTime = new QLineEdit(this);
    m_textDeltaTime = new QLineEdit(this);
    m_textEndTime =   new QLineEdit(this);
    m_scrollBar = new QScrollBar(Qt::Horizontal,this);
    m_scallingBar = new QScrollBar(Qt::Horizontal,this);
    m_simulateButton = new QPushButton("Simulate");
    m_stopButton = new QPushButton("Stop");

    m_showNeButton = new QPushButton("Show Ne");
    m_showEnergyButton = new QPushButton("Show Energy");
    m_showHeavyButton = new QPushButton("Show Heavy Spicies");
    m_showPhiButton = new QPushButton("Show Phi");

    m_vLayoutCheckBoxes->addWidget(m_showNeButton);
    m_vLayoutCheckBoxes->addWidget(m_showEnergyButton);
    m_vLayoutCheckBoxes->addWidget(m_showPhiButton);
    m_vLayoutCheckBoxes->addWidget(m_showHeavyButton);

    m_grid->addWidget(m_customPlot,0,0,2,3);
    m_grid->addWidget(glWidget,2,0,2,3);
    m_grid->addLayout(m_vLayoutCheckBoxes,2,3,2,1);
    m_grid->addWidget(m_simulateButton,6,0);
    m_grid->addWidget(m_stopButton,6,1);
    m_grid->addWidget(m_scallingBar,6,2,1,1);

    m_hLayout->addWidget(new QLabel("Start time:"));
    m_hLayout->addWidget(m_textStartTime);
    m_hLayout->addWidget(new QLabel("dt:"));
    m_hLayout->addWidget(m_textDeltaTime);
    m_hLayout->addWidget(new QLabel("End time:"));
    m_hLayout->addWidget(m_textEndTime);
    m_grid->addLayout(m_hLayout,5,0,1,2);
    m_grid->addWidget(m_scrollBar,5,2,1,1);


    m_data = new simulationData(NX,NY);

    m_data->setDx(1.0/NX);
    m_data->setDy(0.3/NY);
    m_data->setCellsXNumber(NX);
    m_data->setCellsYNumber(NY);
    m_data->setDt(1e-7);

    m_sNe = new solverNe(m_data);
    m_sEn = new solverEnergy(m_data);
    m_sPhi = new solverPhi(m_data);

    m_textStartTime->setText(QString().number(0.0));
    m_textDeltaTime->setText(QString().number(m_data->getDt()));
    m_textEndTime->setText(QString().number(1e-2));

    m_scallingBar->setRange(0,100);
    m_scallingBar->setValue(50);

    connect (m_simulateButton, SIGNAL(clicked(bool)), this, SLOT(simulateData(bool)));
    connect (m_stopButton, SIGNAL(clicked(bool)), this, SLOT(stopAnim(bool)));
    connect(m_scrollBar, SIGNAL(valueChanged(int)), this, SLOT(replotGraph(int)));
    connect(m_scallingBar, SIGNAL(valueChanged(int)), this, SLOT(replotGraph(int)));

    connect(m_showNeButton, SIGNAL(clicked(bool)), this, SLOT(setVisualArrNe()));
    connect(m_showEnergyButton, SIGNAL(clicked(bool)), this, SLOT(setVisualArrEnergy()));
    connect(m_showPhiButton, SIGNAL(clicked(bool)), this, SLOT(setVisualArrPhi()));
    connect(m_showHeavyButton, SIGNAL(clicked(bool)), this, SLOT(setVisualArrHeavy()));




    m_widget->setLayout(m_grid);
    setCentralWidget(m_widget);
    setWindowTitle("PlasmaSolver");
    m_animStopped=true;
    initData();
    //  initData();
}

MainWindow::~MainWindow()
{

}

void MainWindow::replotGraph(int number)
{
    simulationData::simulationParameters* pParams=m_data->getParameters();
    glWidget->setField( m_visualArr , pParams->arrEps, NX, NY, m_data->getDx(), m_data->getDy(),m_scallingBar->value()*m_visualScale);

    glWidget->repaint();//

    /*m_maxY=-1e30;
    m_minY=+1e30;
    m_plots.clear();
    m_plots = m_storage[m_scrollBar->value()].plots;
    m_customPlot->clearGraphs();
    m_customPlot->xAxis->setLabel(QString("x, time= %1").arg(m_time));
    m_customPlot->yAxis->setLabel("y");
    QColor colors[6] = {QColor(255,0,0),QColor(0,255,0),QColor(0,0,255),QColor(255,255,0),QColor(0,255,255),QColor(255,0,255)};
    for (int j = 0; j < m_plots.size(); ++j)
    {
        m_customPlot->addGraph();
        m_customPlot->graph(j)->setName(QString(m_plots[j].name));
        m_customPlot->graph(j)->addToLegend();
        m_customPlot->graph(j)->setData(m_plots[j].x,m_plots[j]. y);
        m_customPlot->graph(j)->setLineStyle((QCPGraph::LineStyle)(1));
        QPen graphPen;
        graphPen.setColor(colors[j]);
        graphPen.setWidthF(2);
        m_customPlot->graph(j)->setPen(graphPen);
        if(!m_checkBoxes[j]->isChecked())
        {
            m_customPlot->graph(j)->setVisible(false);
        }
        else {
            for (int i=0; i < m_plots[j].size; ++i)
            {
                if(m_plots[j].y[i] > m_maxY)
                    m_maxY = m_plots[j].y[i];
                if(m_plots[j].y[i] < m_minY)
                    m_minY = m_plots[j].y[i];
            }

        }
    }
    m_customPlot->xAxis->setRange(-1, 1);
    m_customPlot->yAxis->setRange(m_minY, m_maxY);
    m_customPlot->legend->setVisible(true);
    m_customPlot->replot();*/
}

void MainWindow::saveInStorage()
{
    /*for (int j = 0; j < m_plots.size(); ++j)
    {
        for (int i=0; i < m_plots[j].size; ++i)
        {
            m_plots[j].y[i] = m_plots[j].scale * m_plots[j].arr[i];
        }
    }
    storeStruct storeItem;
    storeItem.plots = m_plots;
    storeItem.time = m_time;
    m_storage.push_back(storeItem);*/

}

void MainWindow::addPlot(double *arr, char *name, int size, double scale)
{
    plotStruct plot;
    plot.x.resize(size);
    plot.y.resize(size);
    plot.arr = arr;
    for (int i = 0; i < size; ++i)
    {
        plot.x[i] = i / (size * 1.0 / 2)-1;
        plot.y[i] = plot.arr[i];
    }
    plot.size = size;
    plot.name = name;
    plot.scale = scale;
    m_plots.push_back(plot);

    QCheckBox* checkBox = new QCheckBox();
    checkBox->setText(QString(name));
    checkBox->setCheckState(Qt::Checked);
    m_checkBoxes.push_back(checkBox);
    m_vLayoutCheckBoxes->addWidget(checkBox);
    //connect(checkBox, SIGNAL(stateChanged(int)), this, SLOT(replotGraph(int)));
}

void MainWindow::addPlotXY(double *arr,double*xx, char *name, int size, double scale)
{
    plotStruct plot;
    plot.x.resize(size);
    plot.y.resize(size);
    plot.arr = arr;
    for (int i = 0; i < size; ++i)
    {
        plot.x[i] = xx[i];//i / (size * 1.0 / 2)-1;
        plot.y[i] = plot.arr[i];
    }
    plot.size = size;
    plot.name = name;
    plot.scale = scale;
    m_plots.push_back(plot);

    QCheckBox* checkBox = new QCheckBox();
    checkBox->setText(QString(name));
    checkBox->setCheckState(Qt::Checked);
    m_checkBoxes.push_back(checkBox);
    m_vLayoutCheckBoxes->addWidget(checkBox);
    //connect(checkBox, SIGNAL(stateChanged(int)), this, SLOT(replotGraph(int)));
}

void MainWindow::initData()
{
    m_fNe = m_data->getFieldNe();
    m_fEnergy = m_data->getFieldEnergy();
    m_fPhi = m_data->getFieldPhi();
    m_numberHeavySpicies = m_data->getNumberHeavySpicies();
    m_time = 0.0;
    for (int j = 0; j < m_numberHeavySpicies; ++j)
    {
        m_fHeavy.push_back(m_data->getFieldHeavySpicies(j));
        m_sHeavy.push_back(new solverHeavySpicies(m_data, j));
    }

    // m_data->setDz(2e-4/NZ);
    //  m_data->setDt(1e-14);

    simulationData::simulationParameters* pParams=m_data->getParameters();
    m_sPhi->init(0.0);
    //m_sNe->init(100000.0);
    for (int i = 0; i < NX; ++i)
    {
        for (int j = 0; j < NY; ++j)
        {
            double x=(i-3*NX/4)*m_data->getDx();
            double y=(j-NY/2)*m_data->getDy();
            double r=x*x+y*y;
            m_fNe->arr[i][j] = 1e5 + 1e11*simulationTools::gauss(sqrt(x*x+y*y), NY*m_data->getDy()*0.1);
            m_fNe->arrPrev[i][j] =m_fNe->arr[i][j];
            m_fEnergy->arr[i][j] = 0.01*m_fNe->arr[i][j];
            m_fEnergy->arrPrev[i][j] = m_fEnergy->arr[i][j];
            for (int h = 0; h < m_numberHeavySpicies; ++h)
            {
                m_fHeavy[h]->arr[i][j] =(m_fNe->arr[i][j]*pParams->T*8.314)/(pParams->p*6.022e23);
                m_fHeavy[h]->arrPrev[i][j] = m_fHeavy[h]->arr[i][j];
            }
        }
    }


    /*for (int i = 0; i < NZ; ++i) {

        double x_=i*m_data->getDz();
        m_fNe->arr[i] =1e5+ 1e11*simulationTools::gauss(x_-1e-4, 2e-5);
        m_fNe->arrPrev[i] =m_fNe->arr[i];
        m_fEnergy->arr[i] = 5.0*m_fNe->arr[i];
        m_fEnergy->arrPrev[i] = m_fEnergy->arr[i];
        m_fPhi->arr[i] = 0.0;
        m_fPhi->arrPrev[i] = 0.0;
        for (int j = 0; j < m_numberHeavySpicies; ++j)
        {
            m_fHeavy[j]->arr[i] =(m_fNe->arr[i]*pParams->T*8.314)/(pParams->p*6.022e23);
            m_fHeavy[j]->arrPrev[i] = m_fHeavy[j]->arr[i];
        }
    }*/
    setVisualArrNe();
    m_data->updateParams();
    replotGraph(0);

    /*m_plots.clear();
    while (QLayoutItem* item = m_vLayoutCheckBoxes->takeAt(0)) {
        delete item->widget();
        delete item;
    }*/
    //addPlot(m_fNe->arr, m_fNe->name ,m_fNe->cellsNumber);
    //addPlot(pParams->arrTe, m_fEnergy->name, m_fEnergy->cellsNumber);
    //addPlot(pParams->arrE, m_fPhi->name, m_fPhi->cellsNumber-1, 1.0);
    for (int j = 0; j < m_numberHeavySpicies; ++j)
    {
        //addPlot(m_fHeavy[j]->arr, m_fHeavy[j]->name, m_fHeavy[j]->cellsNumber);
    }

    /*for (int i = 0; i < m_checkBoxes.size(); ++i) {
        connect(m_checkBoxes[i], SIGNAL(stateChanged(int)), this, SLOT(replotGraph(int)));
    }
    m_scrollBar->setRange(0,1);
    m_storage.clear();
    saveInStorage();*/
    m_sPhi->solve(20);
}

void MainWindow::updateData()
{

    for (int i = 0; i < 1; ++i)
    {
        m_data->updateParams();

       // m_sPhi->solve(2);

        for (int j=0;j<3;j++)
        {
            m_data->calcReaction(simulationData::ReactionName::comsol_eAr_2eArp);
            m_sNe->solve(1);
            m_sEn->solve(1);
            for (int j = 0; j < m_numberHeavySpicies; ++j)
            {
                m_sHeavy[j]->solve(1);
            }
        }


    }
    m_sNe->getStepEuler();
    m_sEn->getStepEuler();
    for (int j = 0; j < m_numberHeavySpicies; ++j)
    {
        m_sHeavy[j]->getStepEuler();
    }

}

void MainWindow::simulateData(bool status)
{
    //if (status==true)
    {

        m_animStopped=false;
        while(!m_animStopped &&  m_time <= (m_textEndTime->text().toDouble() - m_textStartTime->text().toDouble()))
        {
            m_time += 3.0*m_textDeltaTime->text().toDouble();
            for (int i=0;i<3;i++)
            {
                updateData();
            }

            //saveInStorage();
            replotGraph(m_storage.size()!=0 ? m_storage.size()-1 : 0);

            QCoreApplication::processEvents();
        }
    }
}

void MainWindow::setVisualArrNe()
{
    m_visualArr = m_fNe->arr;
    m_visualScale =5e-13;// 0.005;
    replotGraph(m_storage.size()!=0 ? m_storage.size()-1 : 0);
}
void MainWindow::setVisualArrEnergy()
{
    simulationData::simulationParameters* pParams=m_data->getParameters();
    m_visualArr =/*m_data->getReactionRate(simulationData::ReactionName::comsol_eAr_2eArp);*/m_data->getArrTe();// m_fEnergy->arr;
    m_visualScale = 0.01;//2e-12;
    replotGraph(m_storage.size()!=0 ? m_storage.size()-1 : 0);
}
void MainWindow::setVisualArrPhi()
{
    m_visualArr = m_fPhi->arr;
    m_visualScale = 0.00002;
    replotGraph(m_storage.size()!=0 ? m_storage.size()-1 : 0);
}
void MainWindow::setVisualArrHeavy()
{
    m_visualArr = m_fHeavy[0]->arr;
    m_visualScale = 2e-11 * 1e23;
    replotGraph(m_storage.size()!=0 ? m_storage.size()-1 : 0);
}

void MainWindow::stopAnim(bool)
{
    m_animStopped=true;
}



