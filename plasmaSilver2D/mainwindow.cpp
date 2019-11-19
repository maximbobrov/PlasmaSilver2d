
#include "mainwindow.h"
#include "qcustomplot.h"
#include <QGridLayout>
#include <QPushButton>
#include <QtCore>
#include <QVector>
#include <QDebug>
#include <stdio.h>


#define NX 129
#define NY 129



MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    m_widget=new QWidget();

    m_widget->setObjectName("central");
     m_widget->setMinimumSize(900,700);
    glWidget = new GLWidget;



    m_customPlot = new QCustomPlot(this);
    m_grid = new QGridLayout(this);
    m_vLayoutCheckBoxes = new QVBoxLayout(this);
    m_hLayout = new QHBoxLayout(this);
    m_hLayoutButton = new QHBoxLayout(this);
    m_textStartTime = new QLineEdit(this);
    m_textDeltaTime = new QLineEdit(this);
    m_textEndTime =   new QLineEdit(this);
    m_scrollBar = new QScrollBar(Qt::Horizontal,this);
    m_scallingBar = new QScrollBar(Qt::Horizontal,this);
    m_crossBar = new QScrollBar(Qt::Vertical,this);
    m_simulateButton = new QPushButton("Simulate");
    m_stopButton = new QPushButton("Stop");

    m_progressBar = new QProgressBar(this);
    m_progressBar->setRange(0,100);

    m_showNeButton = new QPushButton("Show Ne");
    m_showEnergyButton = new QPushButton("Show Energy");
    m_showHeavyButton = new QPushButton("Show Heavy Spicies");
    m_showPhiButton = new QPushButton("Show Phi");

    m_vLayoutCheckBoxes->addWidget(m_crossBar);
    m_vLayoutCheckBoxes->addWidget(m_showNeButton);
    m_vLayoutCheckBoxes->addWidget(m_showEnergyButton);
    m_vLayoutCheckBoxes->addWidget(m_showPhiButton);
    m_vLayoutCheckBoxes->addWidget(m_showHeavyButton);


    m_grid->addWidget(m_customPlot,0,0,2,3);
    m_grid->addWidget(glWidget,2,0,2,3);
    m_grid->addLayout(m_vLayoutCheckBoxes,2,3,2,1);

    m_hLayout->addWidget(new QLabel("Start time:"));
    m_hLayout->addWidget(m_textStartTime);
    m_hLayout->addWidget(new QLabel("dt:"));
    m_hLayout->addWidget(m_textDeltaTime);
    m_hLayout->addWidget(new QLabel("End time:"));
    m_hLayout->addWidget(m_textEndTime);
    m_grid->addLayout(m_hLayout,5,0,1,2);

    m_hLayoutButton->addWidget(m_simulateButton);
    m_hLayoutButton->addWidget(m_progressBar);
    m_hLayoutButton->addWidget(m_stopButton);
    m_grid->addLayout(m_hLayoutButton,6,0,1,2);

    m_grid->addWidget(m_scrollBar,5,2,1,1);
    m_grid->addWidget(m_scallingBar,6,2,1,1);


    m_data = new simulationData(NX,NY);
    m_visualArr = new double*[NX];
    for (int i = 0; i < NX; ++i)
        m_visualArr[i] = new double [NY];

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
    m_textEndTime->setText(QString().number(1e-4));

    m_scallingBar->setRange(1,100);
    m_scallingBar->setValue(1);

    m_crossBar->setRange(0, NY-1);
    m_crossBar->setValue((NY-1)/2);


    connect (m_simulateButton, SIGNAL(clicked(bool)), this, SLOT(simulateData(bool)));
    connect (m_stopButton, SIGNAL(clicked(bool)), this, SLOT(stopAnim(bool)));
    connect(m_scrollBar, SIGNAL(valueChanged(int)), this, SLOT(replotGraph(int)));
    connect(m_scallingBar, SIGNAL(valueChanged(int)), this, SLOT(setFields(int)));
    connect(m_crossBar, SIGNAL(valueChanged(int)), this, SLOT(setFields(int)));


    connect(m_showNeButton, SIGNAL(clicked(bool)), this, SLOT(setVisualArrNe()));
    connect(m_showEnergyButton, SIGNAL(clicked(bool)), this, SLOT(setVisualArrEnergy()));
    connect(m_showPhiButton, SIGNAL(clicked(bool)), this, SLOT(setVisualArrPhi()));
    connect(m_showHeavyButton, SIGNAL(clicked(bool)), this, SLOT(setVisualArrHeavy()));




    m_widget->setLayout(m_grid);

    setCentralWidget(m_widget);
    setWindowTitle("PlasmaSolver");
    m_animStopped=true;
    initData();
}

MainWindow::~MainWindow()
{

}

void MainWindow::setFields(int number)
{
    replotGraph(m_storage.size()!=0 ? m_storage.size()-1 : 0);
}

void MainWindow::replotGraph(int number)
{
    simulationData::simulationParameters* pParams=m_data->getParameters();

    m_simulateButton->setText("Simulate (time = " + QString().number(m_storage[number].time)+")");

    m_maxY=-1e30;
    m_minY=+1e30;
    m_plots.clear();
    m_plots = m_storage[number].plots;
    m_customPlot->clearGraphs();
    m_customPlot->xAxis->setLabel("x");
    m_customPlot->yAxis->setLabel("y");
    QColor colors[6] = {QColor(255,0,0),QColor(0,255,0),QColor(0,0,255),QColor(255,255,0),QColor(0,255,255),QColor(255,0,255)};
    QVector<double> x, y;
    x.resize(NX);
    y.resize(NX);

    for (int j = 0; j < m_plots.size(); ++j)
    {
        if(m_visualArrName ==  m_plots[j].name)
        {
            for (int ii = 0; ii < NX; ++ii)
                for (int jj = 0; jj < NY; ++jj)
                    m_visualArr[ii][jj] = m_plots[j].arr[ii][jj];

            glWidget->setField( m_visualArr , pParams->arrEps ,pParams->arrMaskPhi, NX, NY, m_data->getDx(), m_data->getDy(), NY-1-m_crossBar->value(), m_scallingBar->value());
            glWidget->repaint();

            m_customPlot->addGraph();
            m_customPlot->graph(0)->setName(m_plots[j].name);
            m_customPlot->graph(0)->addToLegend();
            for (int i = 0; i < NX; ++i)
            {
                x[i] = i / (NX * 1.0 / 2)-1;
                y[i] = m_plots[j].arr[i][(int)(NY-1-m_crossBar->value())];
                if(y[i] > m_maxY)
                    m_maxY = y[i];
                if(y[i] < m_minY)
                    m_minY = y[i];
            }
            m_customPlot->graph(0)->setData(x, y);
            m_customPlot->graph(0)->setLineStyle((QCPGraph::LineStyle)(1));
            QPen graphPen;
            graphPen.setColor(colors[0]);
            graphPen.setWidthF(2);
            m_customPlot->graph(0)->setPen(graphPen);
            /*for (int i=0; i < m_plots[j].size; ++i)
            {
                if(m_plots[j].y[i] > m_maxY)
                    m_maxY = m_plots[j].y[i];
                if(m_plots[j].y[i] < m_minY)
                    m_minY = m_plots[j].y[i];
            }*/

            m_customPlot->xAxis->setRange(-1, 1);
            m_customPlot->yAxis->setRange(m_minY, m_maxY);
            m_customPlot->legend->setVisible(true);
            m_customPlot->replot();
        }


    }
}

void MainWindow::saveInStorage()
{
    for (int k = 0; k < m_plots.size(); ++k)
    {
        for (int i = 0; i < NX; ++i)
        {
            for (int j = 0; j < NY; ++j)
            {
                m_plots[k].arr[i][j] = m_plots[k].arrRef[i][j];
            }
        }
    }
    storeStruct storeItem;
    storeItem.plots = m_plots;
    storeItem.time = m_time;
    m_storage.push_back(storeItem);
}

void MainWindow::addPlot(double **arr, char *name, double scale)
{
    plotStruct plot;
    for (int i = 0; i < NX; ++i)
    {
        for (int j = 0; j < NY; ++j)
        {
            plot.arr[i][j] = arr[i][j];
        }
    }
    plot.arrRef = arr;
    plot.name = name;
    plot.scale = scale;

    m_plots.push_back(plot);
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

    simulationData::simulationParameters* pParams=m_data->getParameters();
    m_sPhi->init(0.0);
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

    m_data->updateParams();
    m_visualArrName = m_fNe->name;

    m_plots.clear();

    addPlot(m_fNe->arr, m_fNe->name);
    addPlot(pParams->arrTe, m_fEnergy->name);
    addPlot(m_fPhi->arr, m_fPhi->name);
    for (int j = 0; j < m_numberHeavySpicies; ++j)
    {
        addPlot(m_fHeavy[j]->arr, m_fHeavy[j]->name);
    }
    m_scrollBar->setRange(0,0);
    m_storage.clear();
    saveInStorage();
    replotGraph(0);
    m_sPhi->solve(20);
}

void MainWindow::updateData(int numberIt)
{

    for (int i = 0; i < numberIt; ++i)
    {
        m_data->updateParams();
        for (int j=0;j<3;j++)
        {
            m_data->calcReaction(simulationData::ReactionName::comsol_eAr_2eArp);
            m_sNe->solve(1);
            m_sEn->solve(1);
            for (int k = 0; k < m_numberHeavySpicies; ++k)
            {
                m_sHeavy[k]->solve(1);
            }
        }

        m_sNe->getStepEuler();
        m_sEn->getStepEuler();
        for (int j = 0; j < m_numberHeavySpicies; ++j)
        {
            m_sHeavy[j]->getStepEuler();
        }
    }

}

void MainWindow::simulateData(bool status)
{
    m_animStopped=false;
    while(!m_animStopped &&  m_time <= (m_textEndTime->text().toDouble() - m_textStartTime->text().toDouble()))
    {
        int saveNum = (m_textEndTime->text().toDouble() - m_textStartTime->text().toDouble()) / m_textDeltaTime->text().toDouble();
        if(saveNum< 20)
            saveNum = 1;
        else
            saveNum = saveNum / 20;
        m_time += m_textStartTime->text().toDouble() + saveNum * m_textDeltaTime->text().toDouble();
        updateData(saveNum);

        saveInStorage();
        m_scrollBar->setRange(0, m_storage.size() - 1);
        m_scrollBar->setValue(m_storage.size() - 1);
        m_simulateButton->setText("Simulate (time = " + QString().number(m_time)+")");
        replotGraph(m_storage.size()!=0 ? m_storage.size()-1 : 0);
        QCoreApplication::processEvents();

        m_progressBar->setValue(100.0 * (m_time - m_textStartTime->text().toDouble()) / (m_textEndTime->text().toDouble() - m_textStartTime->text().toDouble()));

    }
    m_progressBar->setValue(100);
}

void MainWindow::stopAnim(bool)
{
    m_animStopped=true;
}

void MainWindow::setVisualArrNe()
{
    m_visualArrName = m_fNe->name;
    replotGraph(m_storage.size()!=0 ? m_storage.size()-1 : 0);
}
void MainWindow::setVisualArrEnergy()
{
    m_visualArrName = m_fEnergy->name;
    replotGraph(m_storage.size()!=0 ? m_storage.size()-1 : 0);
}
void MainWindow::setVisualArrPhi()
{
    m_visualArrName = m_fPhi->name;
    replotGraph(m_storage.size()!=0 ? m_storage.size()-1 : 0);
}
void MainWindow::setVisualArrHeavy()
{
    m_visualArrName = m_fHeavy[0]->name;
    replotGraph(m_storage.size()!=0 ? m_storage.size()-1 : 0);
}





