#include "column.h"

column::column() {
    //ctor
}

column::~column() {
    //dtor
}

column::column(const column& other) {
    //copy ctor
}

column& column::operator=(const column& rhs) {
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}

column::column(double Dmax) {
    double L,a,aptm;
    L=Dmax;
    if(L<100.0)
    {
        a=0.35*L;
    }
    else
    {
        a=3.48*sqrt(L);
    }
    aptm=0.866025403784439*a;
    vol=3*a*aptm*L;

    CM=Vector3d::Zero();

    N_vertexes=12;
    N_edges=18;
    N_faces=8;

    vertex.resize(3,N_vertexes);

    vertex.col(0)<<-0.5*a,aptm,0.5*L;
    vertex.col(1)<<0.5*a,aptm,0.5*L;
    vertex.col(2)<<a,0.0,0.5*L;
    vertex.col(3)<<0.5*a,-aptm,0.5*L;
    vertex.col(4)<<-0.5*a,-aptm,0.5*L;
    vertex.col(5)<<-a,0.0,0.5*L;
    vertex.col(6)<<-0.5*a,aptm,-0.5*L;
    vertex.col(7)<<0.5*a,aptm,-0.5*L;
    vertex.col(8)<<a,0.0,-0.5*L;
    vertex.col(9)<<0.5*a,-aptm,-0.5*L;
    vertex.col(10)<<-0.5*a,-aptm,-0.5*L;
    vertex.col(11)<<-a,0.0,-0.5*L;

    edge.resize(N_edges);
    for(int i=0;i<N_edges;i++)
    {
        edge[i].resize(2);
    }
    edge[0][0]=0;
    edge[0][1]=1;
    edge[1][0]=1;
    edge[1][1]=2;
    edge[2][0]=2;
    edge[2][1]=3;
    edge[3][0]=3;
    edge[3][1]=4;
    edge[4][0]=4;
    edge[4][1]=5;
    edge[5][0]=5;
    edge[5][1]=0;
    edge[6][0]=6;
    edge[6][1]=7;
    edge[7][0]=7;
    edge[7][1]=8;
    edge[8][0]=8;
    edge[8][1]=9;
    edge[9][0]=9;
    edge[9][1]=10;
    edge[10][0]=10;
    edge[10][1]=11;
    edge[11][0]=11;
    edge[11][1]=6;
    edge[12][0]=0;
    edge[12][1]=6;
    edge[13][0]=1;
    edge[13][1]=7;
    edge[14][0]=2;
    edge[14][1]=8;
    edge[15][0]=3;
    edge[15][1]=9;
    edge[16][0]=4;
    edge[16][1]=10;
    edge[17][0]=5;
    edge[17][1]=11;

    face.resize(N_faces);
    face[0].resize(6);
    face[1].resize(6);
    for(int i=2;i<N_faces;i++)
    {
        face[i].resize(4);
    }
    face[0][0]=5;
    face[0][1]=4;
    face[0][2]=3;
    face[0][3]=2;
    face[0][4]=1;
    face[0][5]=0;

    face[1][0]=6;
    face[1][1]=7;
    face[1][2]=8;
    face[1][3]=9;
    face[1][4]=10;
    face[1][5]=11;

    face[2][0]=0;
    face[2][1]=1;
    face[2][2]=7;
    face[2][3]=6;

    face[3][0]=1;
    face[3][1]=2;
    face[3][2]=8;
    face[3][3]=7;

    face[4][0]=2;
    face[4][1]=3;
    face[4][2]=9;
    face[4][3]=8;

    face[5][0]=3;
    face[5][1]=4;
    face[5][2]=10;
    face[5][3]=9;

    face[6][0]=4;
    face[6][1]=5;
    face[6][2]=11;
    face[6][3]=10;

    face[7][0]=5;
    face[7][1]=0;
    face[7][2]=6;
    face[7][3]=11;

    radius_of_gyration(a/20.);
}

void column::radius_of_gyration(double d)
{
    double Nk,Nj,Ni;
    Ni=-4.*vertex(0,0)/d;
    Nj=2.*vertex(1,0)/d;
    Nk=2.*vertex(2,0)/d;
    double limx=-vertex(0,0);
    double limy1,limy2;
    double offx,offy,offz;
    offx=(floor(Ni)-1)*0.5;
    offy=(floor(Nj)-1)*0.5;
    offz=(floor(Nk)-1)*0.5;
    double x,y,z;
    double R_squared=0;
    int pixels=0;
    for(int k=0;k<Nk;k++)
    {
        z=(k-offz)*d;
        for(int j=0;j<Nj;j++)
        {
            y=(j-offy)*d;
            for(int i=0;i<Ni;i++)
            {
                x=(i-offx)*d;
                limy1=vertex(1,0)*(x/(vertex(0,2))+1.);
                limy2=vertex(1,0)*(x/(vertex(0,2))-1.);
                if(x<-limx)
                {
                    if(abs(y)<limy1)
                    {
                        R_squared+=(x*x+y*y+z*z);
                        pixels++;
                    }
                }
                else if(x>limx)
                {
                    if(abs(y)<limy2)
                    {
                        R_squared+=(x*x+y*y+z*z);
                        pixels++;
                    }

                }
                else if(abs(i)<limx)
                {
                    R_squared+=(x*x+y*y+z*z);
                    pixels++;
                }
            }
        }
    }
    Rgyr=sqrt(R_squared/pixels);
}
