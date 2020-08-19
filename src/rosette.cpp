#include "rosette.h"

rosette6b::rosette6b()
{
    //ctor
}

rosette6b::~rosette6b()
{
    //dtor
}

double rosette6b::f(double x, double k, double l)
{
    return x+k*pow(x,0.63)-l;
};
double rosette6b::df(double x, double k)
{
    return 1+k*pow(x,-0.37);
};

rosette6b::rosette6b(double Dmax)
{
    double len=Dmax*0.5;
    double alpha=28*M_PI/180;
    double tanalpha=tan(alpha);
    double K=(sqrt(3)*1.552)/(2*tanalpha);

/*    inline double f(double x)
    {
        return x+K*pow(x,0.63)-len;
    }
    inline double df(double x)
    {
        return 1+K*0.63*pow(x,-0.37);
    }
    */
    double Lold=0.1;
    double L=0.0;
    double err=1000;
    while(err>0.1)
    {
        //cout<<"NEWTON ";
        //cout<<"  f  "<<f(Lold,K,len)<<"  df  "<<df(Lold,K);
        L=Lold-(f(Lold,K,len))/(df(Lold,K));
        //cout<<"   Lnew "<<L;
        err=abs(L-Lold);
        //cout<<"  err "<<err<<endl;
        Lold=L;
    }
    //cout<<endl<<"L= "<<L<<endl;

    double a=1.552*pow(L,0.63);
    double aptm=0.5*sqrt(3)*a;
    double t=(sqrt(3)*a)/(2*tanalpha);

    vol=3*sqrt(3)*a*a*(3*L+t);
    //cout<<endl<<"vol= "<<vol<<endl;

    CM=Vector3d::Zero();

    N_vertexes=73;
    vertex.resize(3,N_vertexes);


    vertex.col(0)<<0,0,0;

    vertex.col(1)<<-t,a,0;
    vertex.col(2)<<-t,0.5*a,-aptm;
    vertex.col(3)<<-t,-0.5*a,-aptm;
    vertex.col(4)<<-t,-a,0;
    vertex.col(5)<<-t,-0.5*a,aptm;
    vertex.col(6)<<-t,0.5*a,aptm;
    vertex.col(7)<<-len,a,0;
    vertex.col(8)<<-len,0.5*a,-aptm;
    vertex.col(9)<<-len,-0.5*a,-aptm;
    vertex.col(10)<<-len,-a,0;
    vertex.col(11)<<-len,-0.5*a,aptm;
    vertex.col(12)<<-len,0.5*a,aptm;

    vertex.col(13)<<t,a,0;
    vertex.col(14)<<t,0.5*a,-aptm;
    vertex.col(15)<<t,-0.5*a,-aptm;
    vertex.col(16)<<t,-a,0;
    vertex.col(17)<<t,-0.5*a,aptm;
    vertex.col(18)<<t,0.5*a,aptm;
    vertex.col(19)<<len,a,0;
    vertex.col(20)<<len,0.5*a,-aptm;
    vertex.col(21)<<len,-0.5*a,-aptm;
    vertex.col(22)<<len,-a,0;
    vertex.col(23)<<len,-0.5*a,aptm;
    vertex.col(24)<<len,0.5*a,aptm;

    vertex.col(25)<<0,-t,a;
    vertex.col(26)<<-aptm,-t,0.5*a;
    vertex.col(27)<<-aptm,-t,-0.5*a;
    vertex.col(28)<<0,-t,-a;
    vertex.col(29)<<aptm,-t,-0.5*a;
    vertex.col(30)<<aptm,-t,0.5*a;
    vertex.col(31)<<0,-len,a;
    vertex.col(32)<<-aptm,-len,0.5*a;
    vertex.col(33)<<-aptm,-len,-0.5*a;
    vertex.col(34)<<0,-len,-a;
    vertex.col(35)<<aptm,-len,-0.5*a;
    vertex.col(36)<<aptm,-len,0.5*a;

    vertex.col(37)<<0,t,a;
    vertex.col(38)<<-aptm,t,0.5*a;
    vertex.col(39)<<-aptm,t,-0.5*a;
    vertex.col(40)<<0,t,-a;
    vertex.col(41)<<aptm,t,-0.5*a;
    vertex.col(42)<<aptm,t,0.5*a;
    vertex.col(43)<<0,len,a;
    vertex.col(44)<<-aptm,len,0.5*a;
    vertex.col(45)<<-aptm,len,-0.5*a;
    vertex.col(46)<<0,len,-a;
    vertex.col(47)<<aptm,len,-0.5*a;
    vertex.col(48)<<aptm,len,0.5*a;

    vertex.col(49)<<a,0,-t;
    vertex.col(50)<<0.5*a,-aptm,-t;
    vertex.col(51)<<-0.5*a,-aptm,-t;
    vertex.col(52)<<-a,0,-t;
    vertex.col(53)<<-0.5*a,aptm,-t;
    vertex.col(54)<<0.5*a,aptm,-t;
    vertex.col(55)<<a,0,-len;
    vertex.col(56)<<0.5*a,-aptm,-len;
    vertex.col(57)<<-0.5*a,-aptm,-len;
    vertex.col(58)<<-a,0,-len;
    vertex.col(59)<<-0.5*a,aptm,-len;
    vertex.col(60)<<0.5*a,aptm,-len;

    vertex.col(61)<<a,0,t;
    vertex.col(62)<<0.5*a,-aptm,t;
    vertex.col(63)<<-0.5*a,-aptm,t;
    vertex.col(64)<<-a,0,t;
    vertex.col(65)<<-0.5*a,aptm,t;
    vertex.col(66)<<0.5*a,aptm,t;
    vertex.col(67)<<a,0,len;
    vertex.col(68)<<0.5*a,-aptm,len;
    vertex.col(69)<<-0.5*a,-aptm,len;
    vertex.col(70)<<-a,0,len;
    vertex.col(71)<<-0.5*a,aptm,len;
    vertex.col(72)<<0.5*a,aptm,len;

    N_edges=144;
    edge.resize(N_edges);
    for(int i=0;i<N_edges;i++)
    {
        edge[i].resize(2);
    }

    edge[0][0]=0;
    edge[0][1]=1; //Coordinates of bounds of each edge

    edge[1][0]=0;
    edge[1][1]=2;
    edge[2][0]=0;
    edge[2][1]=3;
    edge[3][0]=0;
    edge[3][1]=4;
    edge[4][0]=0;
    edge[4][1]=5;
    edge[5][0]=0;
    edge[5][1]=6;
    edge[6][0]=0;
    edge[6][1]=13;
    edge[7][0]=0;
    edge[7][1]=14;
    edge[8][0]=0;
    edge[8][1]=15;
    edge[9][0]=0;
    edge[9][1]=16;
    edge[10][0]=0;
    edge[10][1]=17;
    edge[11][0]=0;
    edge[11][1]=18;
    
    edge[12][0]=0;
    edge[12][1]=25;
    
    edge[13][0]=0;
    edge[13][1]=26;
    edge[14][0]=0;
    edge[14][1]=27;
    
    edge[15][0]=0;
    edge[15][1]=28;
    
    edge[16][0]=0;
    edge[16][1]=29;
    edge[17][0]=0;
    edge[17][1]=30;
    edge[18][0]=0;
    edge[18][1]=37;
    edge[19][0]=0;
    edge[19][1]=38;
    edge[20][0]=0;
    edge[20][1]=39;
    edge[21][0]=0;
    edge[21][1]=40;
    edge[22][0]=0;
    edge[22][1]=41;
    edge[23][0]=0;
    edge[23][1]=42;
    edge[24][0]=0;
    edge[24][1]=49;
    edge[25][0]=0;
    edge[25][1]=50;
    edge[26][0]=0;
    edge[26][1]=51;
    edge[27][0]=0;
    edge[27][1]=52;
    edge[28][0]=0;
    edge[28][1]=53;
    edge[29][0]=0;
    edge[29][1]=54;
    edge[30][0]=0;
    edge[30][1]=61;
    edge[31][0]=0;
    edge[31][1]=62;
    edge[32][0]=0;
    edge[32][1]=63;
    edge[33][0]=0;
    edge[33][1]=64;
    edge[34][0]=0;
    edge[34][1]=65;
    edge[35][0]=0;
    edge[35][1]=66;

    edge[36][0]=1;
    edge[36][1]=7;
    edge[37][0]=2;
    edge[37][1]=8;
    edge[38][0]=3;
    edge[38][1]=9;
    edge[39][0]=4;
    edge[39][1]=10;
    edge[40][0]=5;
    edge[40][1]=11;
    edge[41][0]=6;
    edge[41][1]=12;
    edge[42][0]=13;
    edge[42][1]=19;
    edge[43][0]=14;
    edge[43][1]=20;
    edge[44][0]=15;
    edge[44][1]=21;
    edge[45][0]=16;
    edge[45][1]=22;
    edge[46][0]=17;
    edge[46][1]=23;
    edge[47][0]=18;
    edge[47][1]=24;
    edge[48][0]=25;
    edge[48][1]=31;
    edge[49][0]=26;
    edge[49][1]=32;
    edge[50][0]=27;
    edge[50][1]=33;
    edge[51][0]=28;
    edge[51][1]=34;
    edge[52][0]=29;
    edge[52][1]=35;
    edge[53][0]=30;
    edge[53][1]=36;
    edge[54][0]=37;
    edge[54][1]=43;
    edge[55][0]=38;
    edge[55][1]=44;
    edge[56][0]=39;
    edge[56][1]=45;
    edge[57][0]=40;
    edge[57][1]=46;
    edge[58][0]=41;
    edge[58][1]=47;
    edge[59][0]=42;
    edge[59][1]=48;
    edge[60][0]=49;
    edge[60][1]=55;
    edge[61][0]=50;
    edge[61][1]=56;
    edge[62][0]=51;
    edge[62][1]=57;
    edge[63][0]=52;
    edge[63][1]=58;
    edge[64][0]=53;
    edge[64][1]=59;
    edge[65][0]=54;
    edge[65][1]=60;
    edge[66][0]=61;
    edge[66][1]=67;
    edge[67][0]=62;
    edge[67][1]=68;
    edge[68][0]=63;
    edge[68][1]=69;
    edge[69][0]=64;
    edge[69][1]=70;
    edge[70][0]=65;
    edge[70][1]=71;
    edge[71][0]=66;
    edge[71][1]=72;
    edge[72][0]=1;
    edge[72][1]=2;
    edge[73][0]=2;
    edge[73][1]=3;
    edge[74][0]=3;
    edge[74][1]=4;
    edge[75][0]=4;
    edge[75][1]=5;
    edge[76][0]=5;
    edge[76][1]=6;
    edge[77][0]=6;
    edge[77][1]=1;
    edge[78][0]=13;
    edge[78][1]=14;
    edge[79][0]=14;
    edge[79][1]=15;
    edge[80][0]=15;
    edge[80][1]=16;
    edge[81][0]=16;
    edge[81][1]=17;
    edge[82][0]=17;
    edge[82][1]=18;
    edge[83][0]=18;
    edge[83][1]=13;
    edge[84][0]=25;
    edge[84][1]=26;
    edge[85][0]=26;
    edge[85][1]=27;
    edge[86][0]=27;
    edge[86][1]=28;
    edge[87][0]=28;
    edge[87][1]=29;
    edge[88][0]=29;
    edge[88][1]=30;
    edge[89][0]=30;
    edge[89][1]=25;
    edge[90][0]=37;
    edge[90][1]=38;
    edge[91][0]=38;
    edge[91][1]=39;
    edge[92][0]=39;
    edge[92][1]=40;
    edge[93][0]=40;
    edge[93][1]=41;
    edge[94][0]=41;
    edge[94][1]=42;
    edge[95][0]=42;
    edge[95][1]=37;
    edge[96][0]=49;
    edge[96][1]=50;
    edge[97][0]=50;
    edge[97][1]=51;
    edge[98][0]=51;
    edge[98][1]=52;
    edge[99][0]=52;
    edge[99][1]=53;
    edge[100][0]=53;
    edge[100][1]=54;
    edge[101][0]=54;
    edge[101][1]=49;
    edge[102][0]=61;
    edge[102][1]=62;
    edge[103][0]=62;
    edge[103][1]=63;
    edge[104][0]=63;
    edge[104][1]=64;
    edge[105][0]=64;
    edge[105][1]=65;
    edge[106][0]=65;
    edge[106][1]=66;
    edge[107][0]=66;
    edge[107][1]=61;
    edge[108][0]=12;
    edge[108][1]=7;
    edge[109][0]=7;
    edge[109][1]=8;
    edge[110][0]=8;
    edge[110][1]=9;
    edge[111][0]=9;
    edge[111][1]=10;
    edge[112][0]=10;
    edge[112][1]=11;
    edge[113][0]=11;
    edge[113][1]=12;
    edge[114][0]=24;
    edge[114][1]=19;
    edge[115][0]=19;
    edge[115][1]=20;
    edge[116][0]=20;
    edge[116][1]=21;
    edge[117][0]=21;
    edge[117][1]=22;
    edge[118][0]=22;
    edge[118][1]=23;
    edge[119][0]=23;
    edge[119][1]=24;
    edge[120][0]=36;
    edge[120][1]=31;
    edge[121][0]=31;
    edge[121][1]=32;
    edge[122][0]=32;
    edge[122][1]=33;
    edge[123][0]=33;
    edge[123][1]=34;
    edge[124][0]=34;
    edge[124][1]=35;
    edge[125][0]=35;
    edge[125][1]=36;
    edge[126][0]=48;
    edge[126][1]=43;
    edge[127][0]=43;
    edge[127][1]=44;
    edge[128][0]=44;
    edge[128][1]=45;
    edge[129][0]=45;
    edge[129][1]=46;
    edge[130][0]=46;
    edge[130][1]=47;
    edge[131][0]=47;
    edge[131][1]=48;
    edge[132][0]=60;
    edge[132][1]=55;
    edge[133][0]=55;
    edge[133][1]=56;
    edge[134][0]=56;
    edge[134][1]=57;
    edge[135][0]=57;
    edge[135][1]=58;
    edge[136][0]=58;
    edge[136][1]=59;
    edge[137][0]=59;
    edge[137][1]=60;
    edge[138][0]=72;
    edge[138][1]=67;
    edge[139][0]=67;
    edge[139][1]=68;
    edge[140][0]=68;
    edge[140][1]=69;
    edge[141][0]=69;
    edge[141][1]=70;
    edge[142][0]=70;
    edge[142][1]=71;
    edge[143][0]=71;
    edge[143][1]=72;

    N_faces=78;
    face.resize(N_faces);
    face[0].resize(6);
    face[1].resize(6);
    face[2].resize(6);
    face[3].resize(6);
    face[4].resize(6);
    face[5].resize(6);
    for(int i=6;i<42;i++)
    {
        face[i].resize(4);
    }
    for(int i=42;i<N_faces;i++)
    {
        face[i].resize(3);
    }

    face[0][0]=7;
    face[0][1]=8;
    face[0][2]=9;
    face[0][3]=10;
    face[0][4]=11;
    face[0][5]=12;

    face[1][0]=72;
    face[1][1]=71;
    face[1][2]=70;
    face[1][3]=69;
    face[1][4]=68;
    face[1][5]=67;

    face[2][0]=48;
    face[2][1]=47;
    face[2][2]=46;
    face[2][3]=45;
    face[2][4]=44;
    face[2][5]=43;

    face[3][0]=55;
    face[3][1]=56;
    face[3][2]=57;
    face[3][3]=58;
    face[3][4]=59;
    face[3][5]=60;

    face[4][0]=24;
    face[4][1]=23;
    face[4][2]=22;
    face[4][3]=21;
    face[4][4]=20;
    face[4][5]=19;

    face[5][0]=31;
    face[5][1]=32;
    face[5][2]=33;
    face[5][3]=34;
    face[5][4]=35;
    face[5][5]=36;

    face[6][0]=24;
    face[6][1]=18;
    face[6][2]=17;
    face[6][3]=23;

    face[7][0]=23;
    face[7][1]=17;
    face[7][2]=16;
    face[7][3]=22;

    face[8][0]=22;
    face[8][1]=16;
    face[8][2]=15;
    face[8][3]=21;

    face[9][0]=21;
    face[9][1]=15;
    face[9][2]=14;
    face[9][3]=20;

    face[10][0]=20;
    face[10][1]=14;
    face[10][2]=13;
    face[10][3]=19;

    face[11][0]=19;
    face[11][1]=13;
    face[11][2]=18;
    face[11][3]=24;

    face[12][0]=6;
    face[12][1]=12;
    face[12][2]=11;
    face[12][3]=5;

    face[13][0]=5;
    face[13][1]=11;
    face[13][2]=10;
    face[13][3]=4;

    face[14][0]=4;
    face[14][1]=10;
    face[14][2]=9;
    face[14][3]=3;

    face[15][0]=3;
    face[15][1]=9;
    face[15][2]=8;
    face[15][3]=2;

    face[16][0]=2;
    face[16][1]=8;
    face[16][2]=7;
    face[16][3]=1;

    face[17][0]=1;
    face[17][1]=7;
    face[17][2]=12;
    face[17][3]=6;

    face[18][0]=48;
    face[18][1]=42;
    face[18][2]=41;
    face[18][3]=47;

    face[19][0]=47;
    face[19][1]=41;
    face[19][2]=40;
    face[19][3]=46;

    face[20][0]=46;
    face[20][1]=40;
    face[20][2]=39;
    face[20][3]=45;

    face[21][0]=45;
    face[21][1]=39;
    face[21][2]=38;
    face[21][3]=44;

    face[22][0]=44;
    face[22][1]=38;
    face[22][2]=37;
    face[22][3]=43;

    face[23][0]=43;
    face[23][1]=37;
    face[23][2]=42;
    face[23][3]=48;

    face[24][0]=36;
    face[24][1]=35;
    face[24][2]=29;
    face[24][3]=30;

    face[25][0]=35;
    face[25][1]=34;
    face[25][2]=28;
    face[25][3]=29;

    face[26][0]=34;
    face[26][1]=33;
    face[26][2]=27;
    face[26][3]=28;

    face[27][0]=33;
    face[27][1]=32;
    face[27][2]=26;
    face[27][3]=27;

    face[28][0]=32;
    face[28][1]=31;
    face[28][2]=25;
    face[28][3]=26;

    face[29][0]=31;
    face[29][1]=36;
    face[29][2]=30;
    face[29][3]=25;

    face[30][0]=72;
    face[30][1]=66;
    face[30][2]=65;
    face[30][3]=71;

    face[31][0]=71;
    face[31][1]=65;
    face[31][2]=64;
    face[31][3]=70;

    face[32][0]=70;
    face[32][1]=64;
    face[32][2]=63;
    face[32][3]=69;

    face[33][0]=69;
    face[33][1]=63;
    face[33][2]=62;
    face[33][3]=68;

    face[34][0]=68;
    face[34][1]=62;
    face[34][2]=61;
    face[34][3]=67;

    face[35][0]=67;
    face[35][1]=61;
    face[35][2]=66;
    face[35][3]=72;

    face[36][0]=54;
    face[36][1]=60;
    face[36][2]=59;
    face[36][3]=53;

    face[37][0]=53;
    face[37][1]=59;
    face[37][2]=58;
    face[37][3]=52;

    face[38][0]=52;
    face[38][1]=58;
    face[38][2]=57;
    face[38][3]=51;

    face[39][0]=51;
    face[39][1]=57;
    face[39][2]=56;
    face[39][3]=50;

    face[40][0]=50;
    face[40][1]=56;
    face[40][2]=55;
    face[40][3]=49;

    face[41][0]=49;
    face[41][1]=55;
    face[41][2]=60;
    face[41][3]=54;


    face[42][0]=26;
    face[42][1]=25;
    face[42][2]=0;

    face[43][0]=25;
    face[43][1]=30;
    face[43][2]=0;

    face[44][0]=30;
    face[44][1]=29;
    face[44][2]=0;

    face[45][0]=29;
    face[45][1]=28;
    face[45][2]=0;

    face[46][0]=28;
    face[46][1]=27;
    face[46][2]=0;

    face[47][0]=27;
    face[47][1]=26;
    face[47][2]=0;

    face[48][0]=37;
    face[48][1]=38;
    face[48][2]=0;

    face[49][0]=38;
    face[49][1]=39;
    face[49][2]=0;

    face[50][0]=39;
    face[50][1]=40;
    face[50][2]=0;

    face[51][0]=40;
    face[51][1]=41;
    face[51][2]=0;

    face[52][0]=41;
    face[52][1]=42;
    face[52][2]=0;

    face[53][0]=42;
    face[53][1]=37;
    face[53][2]=0;

/*
    face[48][0]=37;
    face[48][1]=42;
    face[48][2]=0;

    face[49][0]=38;
    face[49][1]=37;
    face[49][2]=0;

    face[50][0]=39;
    face[50][1]=38;
    face[50][2]=0;

    face[51][0]=40;
    face[51][1]=39;
    face[51][2]=0;

    face[52][0]=41;
    face[52][1]=40;
    face[52][2]=0;

    face[53][0]=42;
    face[53][1]=41;
    face[53][2]=0;
*/
    face[54][0]=13;
    face[54][1]=14;
    face[54][2]=0;

    face[55][0]=14;
    face[55][1]=15;
    face[55][2]=0;

    face[56][0]=15;
    face[56][1]=16;
    face[56][2]=0;

    face[57][0]=16;
    face[57][1]=17;
    face[57][2]=0;

    face[58][0]=17;
    face[58][1]=18;
    face[58][2]=0;

    face[59][0]=18;
    face[59][1]=13;
    face[59][2]=0;

    face[60][0]=2;
    face[60][1]=1;
    face[60][2]=0;

    face[61][0]=1;
    face[61][1]=6;
    face[61][2]=0;

    face[62][0]=6;
    face[62][1]=5;
    face[62][2]=0;

    face[63][0]=5;
    face[63][1]=4;
    face[63][2]=0;

    face[64][0]=4;
    face[64][1]=3;
    face[64][2]=0;

    face[65][0]=3;
    face[65][1]=2;
    face[65][2]=0;

    face[66][0]=50;
    face[66][1]=49;
    face[66][2]=0;

    face[67][0]=49;
    face[67][1]=54;
    face[67][2]=0;

    face[68][0]=54;
    face[68][1]=53;
    face[68][2]=0;

    face[69][0]=53;
    face[69][1]=52;
    face[69][2]=0;

    face[70][0]=52;
    face[70][1]=51;
    face[70][2]=0;

    face[71][0]=51;
    face[71][1]=50;
    face[71][2]=0;
/*
    face[72][0]=61;
    face[72][1]=66;
    face[72][2]=0;

    face[73][0]=62;
    face[73][1]=61;
    face[73][2]=0;

    face[74][0]=63;
    face[74][1]=62;
    face[74][2]=0;

    face[75][0]=64;
    face[75][1]=63;
    face[75][2]=0;

    face[76][0]=65;
    face[76][1]=64;
    face[76][2]=0;

    face[77][0]=66;
    face[77][1]=65;
    face[77][2]=0;
*/
    face[72][0]=61;
    face[72][1]=62;
    face[72][2]=0;

    face[73][0]=62;
    face[73][1]=63;
    face[73][2]=0;

    face[74][0]=63;
    face[74][1]=64;
    face[74][2]=0;

    face[75][0]=64;
    face[75][1]=65;
    face[75][2]=0;

    face[76][0]=65;
    face[76][1]=66;
    face[76][2]=0;

    face[77][0]=66;
    face[77][1]=61;
    face[77][2]=0;
    radius_of_gyration(t/20.);
}

void rosette6b::radius_of_gyration(double d)
{
    MatrixX3i temp;
    fill(temp,d);
    int N_points=temp.rows();
    double R_squared=0.;
    for(int i=0;i<N_points;i++)
    {
        R_squared+=d*d*(temp(i,0)*temp(i,0)+temp(i,1)*temp(i,1)+temp(i,2)*temp(i,2));
    }
    Rgyr=sqrt(R_squared/N_points);
}
