#include "dendrite.h"

dendrite::dendrite()
{
    //ctor
}

dendrite::~dendrite()
{
    //dtor
}

dendrite::dendrite(double Dmax)
{
    double a=0.5*Dmax;
    double L;
    if(a<=2.0)
    {
        L=Dmax;
    }
    else if(a>=5.0)
    {
        L=2.4883*pow(a,0.474);
    }
    else
    {
        L=2+0.25*(2.4883*pow(a,0.474)-2)*(a-1);
    }

    CM=Vector3d::Zero();

    N_vertexes=72;

    vertex.resize(3,N_vertexes);

    double len=Dmax*0.5;
    double len_med=len*0.3;
    double len_bra=len*0.8;
    double len_min=len_med*0.6;
    double h=L;
    double aptm=sqrt(3)*len*0.5;
    double aptm_min=sqrt(3)*len_min*0.5;
    double aptm_bra=sqrt(3)*len_bra*0.5;
    double theta=M_PI/3.;

    vol=3.*len_min*aptm_min*h +6.*(
          h*len_min*(len_med+0.5*len_min)
         +2.*(h*len_min*0.5*(len_bra-len_min))/3.
         +h*len_min*(len-(aptm_min+len_med+0.5*len_min))/3.
                                  );

    vertex.col(0)<<0,len,0;
    vertex.col(1)<<aptm,0.5*len,0;
    vertex.col(2)<<aptm,-0.5*len,0;
    vertex.col(3)<<0,-len,0;
    vertex.col(4)<<-aptm,-0.5*len,0;
    vertex.col(5)<<-aptm,0.5*len,0;

    vertex.col(6)<<-0.5*len_bra,aptm_bra,0;
    vertex.col(7)<<0.5*len_bra,aptm_bra,0;
    vertex.col(8)<<len_bra,0,0;
    vertex.col(9)<<0.5*len_bra,-aptm_bra,0;
    vertex.col(10)<<-0.5*len_bra,-aptm_bra,0;
    vertex.col(11)<<-len_bra,0,0;

    vertex.col(12)<<-len_min*0.5,aptm_min,-0.5*h;
    vertex.col(13)<<len_min*0.5,aptm_min,-0.5*h;
    vertex.col(14)<<len_min,0,-0.5*h;
    vertex.col(15)<<len_min*0.5,-aptm_min,-0.5*h;
    vertex.col(16)<<-len_min*0.5,-aptm_min,-0.5*h;
    vertex.col(17)<<-len_min,0,-0.5*h;
    vertex.col(18)<<-len_min*0.5,aptm_min,0.5*h;
    vertex.col(19)<<len_min*0.5,aptm_min,0.5*h;
    vertex.col(20)<<len_min,0,0.5*h;
    vertex.col(21)<<len_min*0.5,-aptm_min,0.5*h;
    vertex.col(22)<<-len_min*0.5,-aptm_min,0.5*h;
    vertex.col(23)<<-len_min,0,0.5*h;

    vertex.col(24)<<vertex(0,12),vertex(1,12)+len_med+0.5*len_min,vertex(2,12);
    vertex.col(25)<<vertex(0,12),vertex(1,12)+len_med-0.5*len_min,vertex(2,12);
    vertex.col(26)<<vertex(0,13),vertex(1,13)+len_med+0.5*len_min,vertex(2,13);
    vertex.col(27)<<vertex(0,13),vertex(1,13)+len_med-0.5*len_min,vertex(2,13);
    vertex.col(28)<<vertex(0,24)*cos(theta)+vertex(1,24)*sin(theta),-vertex(0,24)*sin(theta)+vertex(1,24)*cos(theta),-0.5*h;
    vertex.col(29)<<vertex(0,25)*cos(theta)+vertex(1,25)*sin(theta),-vertex(0,25)*sin(theta)+vertex(1,25)*cos(theta),-0.5*h;
    vertex.col(30)<<vertex(0,26)*cos(theta)+vertex(1,26)*sin(theta),-vertex(0,26)*sin(theta)+vertex(1,26)*cos(theta),-0.5*h;
    vertex.col(31)<<vertex(0,27)*cos(theta)+vertex(1,27)*sin(theta),-vertex(0,27)*sin(theta)+vertex(1,27)*cos(theta),-0.5*h;
    vertex.col(32)<<vertex(0,24)*cos(theta)-vertex(1,24)*sin(theta),vertex(0,24)*sin(theta)+vertex(1,24)*cos(theta),-0.5*h;
    vertex.col(33)<<vertex(0,25)*cos(theta)-vertex(1,25)*sin(theta),vertex(0,25)*sin(theta)+vertex(1,25)*cos(theta),-0.5*h;
    vertex.col(34)<<vertex(0,26)*cos(theta)-vertex(1,26)*sin(theta),vertex(0,26)*sin(theta)+vertex(1,26)*cos(theta),-0.5*h;
    vertex.col(35)<<vertex(0,27)*cos(theta)-vertex(1,27)*sin(theta),vertex(0,27)*sin(theta)+vertex(1,27)*cos(theta),-0.5*h;
    vertex.col(36)<<vertex(0,24),-vertex(1,24),-0.5*h;
    vertex.col(37)<<vertex(0,25),-vertex(1,25),-0.5*h;
    vertex.col(38)<<vertex(0,26),-vertex(1,26),-0.5*h;
    vertex.col(39)<<vertex(0,27),-vertex(1,27),-0.5*h;
    vertex.col(40)<<vertex(0,28),-vertex(1,28),-0.5*h;
    vertex.col(41)<<vertex(0,29),-vertex(1,29),-0.5*h;
    vertex.col(42)<<vertex(0,30),-vertex(1,30),-0.5*h;
    vertex.col(43)<<vertex(0,31),-vertex(1,31),-0.5*h;
    vertex.col(44)<<vertex(0,32),-vertex(1,32),-0.5*h;
    vertex.col(45)<<vertex(0,33),-vertex(1,33),-0.5*h;
    vertex.col(46)<<vertex(0,34),-vertex(1,34),-0.5*h;
    vertex.col(47)<<vertex(0,35),-vertex(1,35),-0.5*h;

    vertex.col(48)<<vertex(0,24),vertex(1,24),0.5*h;
    vertex.col(49)<<vertex(0,25),vertex(1,25),0.5*h;
    vertex.col(50)<<vertex(0,26),vertex(1,26),0.5*h;
    vertex.col(51)<<vertex(0,27),vertex(1,27),0.5*h;
    vertex.col(52)<<vertex(0,28),vertex(1,28),0.5*h;
    vertex.col(53)<<vertex(0,29),vertex(1,29),0.5*h;
    vertex.col(54)<<vertex(0,30),vertex(1,30),0.5*h;
    vertex.col(55)<<vertex(0,31),vertex(1,31),0.5*h;
    vertex.col(56)<<vertex(0,32),vertex(1,32),0.5*h;
    vertex.col(57)<<vertex(0,33),vertex(1,33),0.5*h;
    vertex.col(58)<<vertex(0,34),vertex(1,34),0.5*h;
    vertex.col(59)<<vertex(0,35),vertex(1,35),0.5*h;
    vertex.col(60)<<vertex(0,36),vertex(1,36),0.5*h;
    vertex.col(61)<<vertex(0,37),vertex(1,37),0.5*h;
    vertex.col(62)<<vertex(0,38),vertex(1,38),0.5*h;
    vertex.col(63)<<vertex(0,39),vertex(1,39),0.5*h;
    vertex.col(64)<<vertex(0,40),vertex(1,40),0.5*h;
    vertex.col(65)<<vertex(0,41),vertex(1,41),0.5*h;
    vertex.col(66)<<vertex(0,42),vertex(1,42),0.5*h;
    vertex.col(67)<<vertex(0,43),vertex(1,43),0.5*h;
    vertex.col(68)<<vertex(0,44),vertex(1,44),0.5*h;
    vertex.col(69)<<vertex(0,45),vertex(1,45),0.5*h;
    vertex.col(70)<<vertex(0,46),vertex(1,46),0.5*h;
    vertex.col(71)<<vertex(0,47),vertex(1,47),0.5*h;

    N_edges=96;
    edge.resize(N_edges);
    for(int i=0;i<N_edges;i++)
    {
        edge[i].resize(2);
    }
    edge[0][0]=12;
    edge[0][1]=24;
    edge[1][0]=13;
    edge[1][1]=26;
    edge[2][0]=13;
    edge[2][1]=28;
    edge[3][0]=14;
    edge[3][1]=30;
    edge[4][0]=14;
    edge[4][1]=42;
    edge[5][0]=15;
    edge[5][1]=40;
    edge[6][0]=15;
    edge[6][1]=38;
    edge[7][0]=16;
    edge[7][1]=36;
    edge[8][0]=16;
    edge[8][1]=46;
    edge[9][0]=17;
    edge[9][1]=44;
    edge[10][0]=17;
    edge[10][1]=32;
    edge[11][0]=12;
    edge[11][1]=34;
    edge[12][0]=18;
    edge[12][1]=48;
    edge[13][0]=19;
    edge[13][1]=50;
    edge[14][0]=19;
    edge[14][1]=52;
    edge[15][0]=20;
    edge[15][1]=54;
    edge[16][0]=20;
    edge[16][1]=66;
    edge[17][0]=21;
    edge[17][1]=64;
    edge[18][0]=21;
    edge[18][1]=62;
    edge[19][0]=22;
    edge[19][1]=60;
    edge[20][0]=22;
    edge[20][1]=70;
    edge[21][0]=23;
    edge[21][1]=68;
    edge[22][0]=23;
    edge[22][1]=56;
    edge[23][0]=18;
    edge[23][1]=58;
    edge[24][0]=24;
    edge[24][1]=0;
    edge[25][0]=48;
    edge[25][1]=0;
    edge[26][0]=26;
    edge[26][1]=0;
    edge[27][0]=50;
    edge[27][1]=0;
    edge[28][0]=28;
    edge[28][1]=1;
    edge[29][0]=30;
    edge[29][1]=1;
    edge[30][0]=52;
    edge[30][1]=1;
    edge[31][0]=54;
    edge[31][1]=1;
    edge[32][0]=40;
    edge[32][1]=2;
    edge[33][0]=42;
    edge[33][1]=2;
    edge[34][0]=64;
    edge[34][1]=2;
    edge[35][0]=66;
    edge[35][1]=2;
    edge[36][0]=36;
    edge[36][1]=3;
    edge[37][0]=38;
    edge[37][1]=3;
    edge[38][0]=60;
    edge[38][1]=3;
    edge[39][0]=62;
    edge[39][1]=3;
    edge[40][0]=44;
    edge[40][1]=4;
    edge[41][0]=46;
    edge[41][1]=4;
    edge[42][0]=70;
    edge[42][1]=4;
    edge[43][0]=68;
    edge[43][1]=4;
    edge[44][0]=32;
    edge[44][1]=5;
    edge[45][0]=34;
    edge[45][1]=5;
    edge[46][0]=56;
    edge[46][1]=5;
    edge[47][0]=58;
    edge[47][1]=5;
    edge[48][0]=34;
    edge[48][1]=6;
    edge[49][0]=58;
    edge[49][1]=6;
    edge[50][0]=35;
    edge[50][1]=6;
    edge[51][0]=59;
    edge[51][1]=6;
    edge[52][0]=25;
    edge[52][1]=6;
    edge[53][0]=49;
    edge[53][1]=6;
    edge[54][0]=24;
    edge[54][1]=6;
    edge[55][0]=48;
    edge[55][1]=6;
    edge[56][0]=27;
    edge[56][1]=7;
    edge[57][0]=51;
    edge[57][1]=7;
    edge[58][0]=26;
    edge[58][1]=7;
    edge[59][0]=50;
    edge[59][1]=7;
    edge[60][0]=29;
    edge[60][1]=7;
    edge[61][0]=53;
    edge[61][1]=7;
    edge[62][0]=28;
    edge[62][1]=7;
    edge[63][0]=52;
    edge[63][1]=7;
    edge[64][0]=31;
    edge[64][1]=8;
    edge[65][0]=55;
    edge[65][1]=8;
    edge[66][0]=30;
    edge[66][1]=8;
    edge[67][0]=54;
    edge[67][1]=8;
    edge[68][0]=43;
    edge[68][1]=8;
    edge[69][0]=67;
    edge[69][1]=8;
    edge[70][0]=42;
    edge[70][1]=8;
    edge[71][0]=66;
    edge[71][1]=8;
    edge[72][0]=65;
    edge[72][1]=9;
    edge[73][0]=41;
    edge[73][1]=9;
    edge[74][0]=40;
    edge[74][1]=9;
    edge[75][0]=64;
    edge[75][1]=9;
    edge[76][0]=63;
    edge[76][1]=9;
    edge[77][0]=39;
    edge[77][1]=9;
    edge[78][0]=38;
    edge[78][1]=9;
    edge[79][0]=62;
    edge[79][1]=9;
    edge[80][0]=61;
    edge[80][1]=10;
    edge[81][0]=37;
    edge[81][1]=10;
    edge[82][0]=60;
    edge[82][1]=10;
    edge[83][0]=36;
    edge[83][1]=10;
    edge[84][0]=47;
    edge[84][1]=10;
    edge[85][0]=71;
    edge[85][1]=10;
    edge[86][0]=46;
    edge[86][1]=10;
    edge[87][0]=70;
    edge[87][1]=10;
    edge[88][0]=69;
    edge[88][1]=11;
    edge[89][0]=45;
    edge[89][1]=11;
    edge[90][0]=44;
    edge[90][1]=11;
    edge[91][0]=68;
    edge[91][1]=11;
    edge[92][0]=33;
    edge[92][1]=11;
    edge[93][0]=57;
    edge[93][1]=11;
    edge[94][0]=32;
    edge[94][1]=11;
    edge[95][0]=56;
    edge[95][1]=11;

    N_faces=98;
    face.resize(N_faces);
    face[0].resize(6);
    face[1].resize(6);
    for(int i=2;i<26;i++)
    {
        face[i].resize(4);
    }
    for(int i=26;i<N_faces;i++)
    {
        face[i].resize(3);
    }

    face[0][0]=23;
    face[0][1]=22;
    face[0][2]=21;
    face[0][3]=20;
    face[0][4]=19;
    face[0][5]=18;

    face[1][0]=12;
    face[1][1]=13;
    face[1][2]=14;
    face[1][3]=15;
    face[1][4]=16;
    face[1][5]=17;

    face[2][0]=12;
    face[2][1]=18;
    face[2][2]=49;
    face[2][3]=25;

    face[3][0]=18;
    face[3][1]=19;
    face[3][2]=50;
    face[3][3]=48;

    face[4][0]=24;
    face[4][1]=26;
    face[4][2]=13;
    face[4][3]=12;

    face[5][0]=27;
    face[5][1]=51;
    face[5][2]=19;
    face[5][3]=13;

    face[6][0]=23;
    face[6][1]=18;
    face[6][2]=58;
    face[6][3]=56;

    face[7][0]=32;
    face[7][1]=34;
    face[7][2]=12;
    face[7][3]=17;

    face[8][0]=17;
    face[8][1]=23;
    face[8][2]=57;
    face[8][3]=33;

    face[9][0]=35;
    face[9][1]=59;
    face[9][2]=18;
    face[9][3]=12;

    face[10][0]=22;
    face[10][1]=23;
    face[10][2]=68;
    face[10][3]=70;

    face[11][0]=46;
    face[11][1]=44;
    face[11][2]=17;
    face[11][3]=16;

    face[12][0]=16;
    face[12][1]=22;
    face[12][2]=71;
    face[12][3]=47;

    face[13][0]=45;
    face[13][1]=69;
    face[13][2]=23;
    face[13][3]=17;

    face[14][0]=21;
    face[14][1]=22;
    face[14][2]=60;
    face[14][3]=62;

    face[15][0]=38;
    face[15][1]=36;
    face[15][2]=16;
    face[15][3]=15;

    face[16][0]=15;
    face[16][1]=21;
    face[16][2]=63;
    face[16][3]=39;

    face[17][0]=61;
    face[17][1]=22;
    face[17][2]=16;
    face[17][3]=37;

    face[18][0]=20;
    face[18][1]=21;
    face[18][2]=64;
    face[18][3]=66;

    face[19][0]=40;
    face[19][1]=15;
    face[19][2]=14;
    face[19][3]=42;

    face[20][0]=20;
    face[20][1]=67;
    face[20][2]=43;
    face[20][3]=14;

    face[21][0]=65;
    face[21][1]=21;
    face[21][2]=15;
    face[21][3]=41;

    face[22][0]=19;
    face[22][1]=20;
    face[22][2]=54;
    face[22][3]=52;

    face[23][0]=28;
    face[23][1]=30;
    face[23][2]=14;
    face[23][3]=13;

    face[24][0]=13;
    face[24][1]=19;
    face[24][2]=53;
    face[24][3]=29;

    face[25][0]=55;
    face[25][1]=20;
    face[25][2]=14;
    face[25][3]=31;

    face[26][0]=48;
    face[26][1]=50;
    face[26][2]=0;

    face[27][0]=24;
    face[27][1]=0;
    face[27][2]=26;

    face[28][0]=48;
    face[28][1]=0;
    face[28][2]=24;

    face[29][0]=0;
    face[29][1]=50;
    face[29][2]=26;

    face[30][0]=6;
    face[30][1]=48;
    face[30][2]=24;

    face[31][0]=6;
    face[31][1]=49;
    face[31][2]=48;

    face[32][0]=6;
    face[32][1]=25;
    face[32][2]=49;

    face[33][0]=6;
    face[33][1]=24;
    face[33][2]=25;

    face[34][0]=7;
    face[34][1]=50;
    face[34][2]=51;

    face[35][0]=7;
    face[35][1]=26;
    face[35][2]=50;

    face[36][0]=7;
    face[36][1]=51;
    face[36][2]=27;

    face[37][0]=7;
    face[37][1]=27;
    face[37][2]=26;

    face[38][0]=1;
    face[38][1]=52;
    face[38][2]=54;

    face[39][0]=1;
    face[39][1]=28;
    face[39][2]=52;

    face[40][0]=1;
    face[40][1]=30;
    face[40][2]=28;

    face[41][0]=1;
    face[41][1]=54;
    face[41][2]=30;

    face[42][0]=7;
    face[42][1]=53;
    face[42][2]=52;

    face[43][0]=7;
    face[43][1]=52;
    face[43][2]=28;

    face[44][0]=7;
    face[44][1]=28;
    face[44][2]=29;

    face[45][0]=7;
    face[45][1]=29;
    face[45][2]=53;

    face[46][0]=8;
    face[46][1]=54;
    face[46][2]=55;

    face[47][0]=8;
    face[47][1]=55;
    face[47][2]=31;

    face[48][0]=8;
    face[48][1]=31;
    face[48][2]=30;

    face[49][0]=8;
    face[49][1]=30;
    face[49][2]=54;

    face[50][0]=2;
    face[50][1]=66;
    face[50][2]=64;

    face[51][0]=2;
    face[51][1]=42;
    face[51][2]=66;

    face[52][0]=2;
    face[52][1]=64;
    face[52][2]=40;

    face[53][0]=2;
    face[53][1]=40;
    face[53][2]=42;

    face[54][0]=8;
    face[54][1]=67;
    face[54][2]=66;

    face[55][0]=8;
    face[55][1]=66;
    face[55][2]=42;

    face[56][0]=8;
    face[56][1]=43;
    face[56][2]=67;

    face[57][0]=8;
    face[57][1]=42;
    face[57][2]=43;

    face[58][0]=9;
    face[58][1]=64;
    face[58][2]=65;

    face[59][0]=9;
    face[59][1]=40;
    face[59][2]=64;

    face[60][0]=9;
    face[60][1]=41;
    face[60][2]=40;

    face[61][0]=9;
    face[61][1]=65;
    face[61][2]=41;

    face[62][0]=3;
    face[62][1]=62;
    face[62][2]=60;

    face[63][0]=3;
    face[63][1]=38;
    face[63][2]=62;

    face[64][0]=3;
    face[64][1]=60;
    face[64][2]=36;

    face[65][0]=3;
    face[65][1]=36;
    face[65][2]=38;

    face[66][0]=9;
    face[66][1]=63;
    face[66][2]=62;

    face[67][0]=9;
    face[67][1]=62;
    face[67][2]=38;

    face[68][0]=9;
    face[68][1]=38;
    face[68][2]=39;

    face[69][0]=9;
    face[69][1]=39;
    face[69][2]=63;

    face[70][0]=10;
    face[70][1]=60;
    face[70][2]=61;

    face[71][0]=10;
    face[71][1]=61;
    face[71][2]=37;

    face[72][0]=10;
    face[72][1]=37;
    face[72][2]=36;

    face[73][0]=10;
    face[73][1]=36;
    face[73][2]=60;

    face[74][0]=4;
    face[74][1]=70;
    face[74][2]=68;

    face[75][0]=4;
    face[75][1]=46;
    face[75][2]=70;

    face[76][0]=4;
    face[76][1]=44;
    face[76][2]=46;

    face[77][0]=4;
    face[77][1]=68;
    face[77][2]=44;

    face[78][0]=10;
    face[78][1]=71;
    face[78][2]=70;

    face[79][0]=10;
    face[79][1]=70;
    face[79][2]=46;

    face[80][0]=10;
    face[80][1]=46;
    face[80][2]=47;

    face[81][0]=10;
    face[81][1]=47;
    face[81][2]=71;

    face[82][0]=11;
    face[82][1]=68;
    face[82][2]=69;

    face[83][0]=11;
    face[83][1]=69;
    face[83][2]=45;

    face[84][0]=11;
    face[84][1]=45;
    face[84][2]=44;

    face[85][0]=11;
    face[85][1]=44;
    face[85][2]=68;

    face[86][0]=5;
    face[86][1]=56;
    face[86][2]=58;

    face[87][0]=5;
    face[87][1]=58;
    face[87][2]=34;

    face[88][0]=5;
    face[88][1]=32;
    face[88][2]=56;

    face[89][0]=5;
    face[89][1]=34;
    face[89][2]=32;

    face[90][0]=11;
    face[90][1]=57;
    face[90][2]=56;

    face[91][0]=11;
    face[91][1]=56;
    face[91][2]=32;

    face[92][0]=11;
    face[92][1]=32;
    face[92][2]=33;

    face[93][0]=11;
    face[93][1]=33;
    face[93][2]=57;

    face[94][0]=6;
    face[94][1]=58;
    face[94][2]=59;

    face[95][0]=6;
    face[95][1]=59;
    face[95][2]=35;

    face[96][0]=6;
    face[96][1]=35;
    face[96][2]=34;

    face[97][0]=6;
    face[97][1]=34;
    face[97][2]=58;

    radius_of_gyration(h/20.);
}

void dendrite::radius_of_gyration(double d)
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
