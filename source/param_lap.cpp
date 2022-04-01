/*平面パラメータ化
  実行までの例（Eigenのパスは自身のパスを指定してください）
  g++ param_lap.cpp -o param_lap -std=c++17 -I /usr/local/.../include/eigen3
  ./param_lap sample.off
*/
# include<iostream>
# include<cmath>
# include<cstdio>
# include<vector>
# include<algorithm>
# include<Eigen/Eigen>
# include<fstream>
# define PI 3.14159265358979323846
# define N 100
/*using SparseMatrix = Eigen::SparseMatrix<double>;
using SparseSolver = Eigen::SparseLU<SparseMatrix>;
using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;*/
using namespace std;
using namespace Eigen;
struct edge{
  int v1;
  int v2;
  int idx;
};
struct face{
  int v1;
  int v2;
  int v3;
};
struct point{
  double x;
  double y;
  double z;
};
edge assain(int a, int b, int c){
  if(a < b) return {a, b, c};
  else return {b, a, c};
}
int calc_v_tgt(face f, int half){
  int i = ((half % 3) + 2) % 3;
  if(i == 0) return f.v1;
  else if(i == 1) return f.v2;
  else return f.v3;
  return 0;
}
int calc_opposite_v(face f, int half){
  int i = half % 3;
  if(i == 0) return f.v1;
  else if(i == 1) return f.v2;
  else return f.v3;
  return 0;
}
//cot重みを計算する関数
double calc_cot_weight(point v1, point v2, point v3, point v4){
  double cos_v1, cos_v2, sin_v1, sin_v2;
  Vector3d a, b, c, d;
  a << v1.x-v2.x, v1.y-v2.y, v1.z-v2.z;
  b << v3.x-v2.x, v3.y-v2.y, v3.z-v2.z;
  c << v3.x-v4.x, v3.y-v4.y, v3.z-v4.z;
  d << v1.x-v4.x, v1.y-v4.y, v1.z-v4.z;
  
  double l_ab, l_cd;
  l_ab = a.norm() * b.norm();
  l_cd = c.norm() * d.norm();

  cos_v1 = a.dot(b) / l_ab;
  cos_v2 = c.dot(d) / l_cd;
  if(cos_v1 < 0) cout <<"cos_v1 error: " << cos_v1 <<"\n";
  if(cos_v2 < 0) cout <<"cos_v2 error: " << cos_v2 <<"\n";

  sin_v1 = ((a.cross(b)).norm()) / l_ab;
  sin_v2 = ((c.cross(d)).norm()) / l_cd;
  if(sin_v1 < 0) cout <<"sin_v1 error: " << sin_v1 <<"\n";
  if(sin_v2 < 0) cout <<"sin_v2 error: " << sin_v2 <<"\n";

  double w = cos_v1/sin_v1 + cos_v2/sin_v2;

  return w;
}

int main(int arg, char *fname[]){
  FILE* fl = fopen(fname[1], "r");
  if(!fl){
    perror("File opening error");
    return 0;
  }

//halfedge  /////////////////////////////////////////////////////////////////
  printf("halfedge\n");
  int v_num, f_num, t;
  for(int i=0; i<2; i++){
    char c[N] = {0};
    fgets(c, N, fl);
    if(i == 1) sscanf(c, "%d %d %d", &v_num, &f_num, &t);
  }

  vector<double> z;
  vector<point> p;
  for(int i=0; i<v_num; i++){
    char c[N] = {0};
    fgets(c, N, fl);
    double x, y, zz;
    sscanf(c, "%lf %lf %lf", &x, &y, &zz);
    p.push_back({x,y,zz});
    z.push_back(zz);
  }

  vector<edge> A;
  vector<face> F;
  vector<int> outgoing_halfedges(v_num, -1);
  for(int i=0; i<f_num; i++){
    char c[N] = {0};
    fgets(c, N, fl);
    int a[3];
    sscanf(c, "%d %d %d %d", &t, &a[0], &a[1], &a[2]);
    face f = {a[0], a[1], a[2]};
    F.push_back(f);
    for(int j=0; j<3; j++){
      edge e = assain(a[j], a[(j+1)%3], 3*i+((j+2)%3));
      A.push_back(e);
      if(outgoing_halfedges[a[j]] == -1) outgoing_halfedges[a[j]] = 3*i+((j+2)%3);
    }
  }
/*  printf("face\n");
  for(auto out : F){
    printf("%4d        %4d        %4d\n",out.v1, out.v2, out.v3);
  }*/

/*  printf("outgoing_halfedges\nvertex_idx  edge_idx\n");
  int cnt = 0;
  for(auto out : outgoing_halfedges){
    printf("%4d        %4d\n",cnt, out);
    cnt++;
  }*/

  sort(A.begin(), A.end(),
      [](const edge& x, const edge& y){
      if(x.v1 == y.v1) return x.v2 < y.v2;
      else return x.v1 < y.v1;});

/*  printf("\nsort\n");
  for(auto out : A){
    printf("%d  %d  %d\n", out.v1, out.v2, out.idx);
  };*/

  vector<int> opposite_halfedges(f_num*3, -1);
  int ii=0;
  while(ii < f_num*3){
    if(A[ii].v1 == A[ii+1].v1 && A[ii].v2 == A[ii+1].v2){
      opposite_halfedges[A[ii].idx] = A[ii+1].idx;
      opposite_halfedges[A[ii+1].idx] = A[ii].idx;
      ii+=2;
    }else ii++;
  }
/*  printf("\nopposite_halfedges\nopposite_idx  edge_idx\n");
  for(auto out : opposite_halfedges){
    if(out > -1) printf("%4d         %4d\n", out, opposite_halfedges[out]);
    else printf("%4d\n", out);
  };
  cout << v_num <<"\n";*/
//メッシュのパラメータ化 //////////////////////////////////////////////////////////
  printf("メッシュのパラメータ化\n");
  MatrixXd laplacian(v_num, v_num); laplacian = MatrixXd::Zero(v_num, v_num);
  VectorXd bx(v_num); bx = VectorXd::Zero(v_num);
  VectorXd by(v_num); by = VectorXd::Zero(v_num);

/*  SparseMatrix laplacian(v_num, v_num);
  Vector bx = Vector::Zero(v_num);
  Vector by = Vector::Zero(v_num);*/

  int v_baundray = 0;
  int bundray_cnt = 0;
  for(int i=0; i<v_num; i++){ //点iに接続している頂点数の算出し、laplacianを求める
    int h = outgoing_halfedges[i]; //点iを支点とした半辺の添字
    int h_end = h; //終了条件（半時計周りに半辺を通して頂点を調べるためスタート地点を保存する）
    //頂点数をで割って重みを求める部分をcot重みに変更する。
    double sum = 0; //繋がっている頂点の重みの合計
    do{
      int j = calc_v_tgt(F[h/3], h);
      //-1をcot重みに変更する。４点渡して求める。calc_cot_weight関数を利用。
      point v1, v2, v3, v4;
      v1 = p[i];
      v2 = p[calc_opposite_v(F[h/3], h)];
      v3 = p[j];
      int ophf = opposite_halfedges[h];
      v4 = p[calc_opposite_v(F[ophf/3], ophf)];

      double w = calc_cot_weight(v1,v2,v3,v4);
      if (w < 0) cout <<"error w: " << w <<"\n";
      laplacian(i, j) = w;
      sum += w;
      int h_prev = 0;
      if(h%3 != 0) h_prev = h-1;
      else h_prev = h+2;
      h = opposite_halfedges[h_prev];
      if(h == -1){ //境界辺に当たった場合の処理
        sum = -1;
        bundray_cnt++;
        v_baundray = i;
        for(int k=0; k<v_num; k++){
          laplacian(i, k) = 0;
        }
        laplacian(i, i) = 1;
        break;
      }
    }while(h != h_end);
    if(sum != -1) laplacian(i, i) = -1*sum;
  }
// 境界面の処理（初期値となる多角形の外周の値をbx、byに代入）
  printf("境界面の処理\n");
  //cout << bundray_cnt <<"\n";
  vector<point> polygon;
  for(int i=0; i<bundray_cnt; i++){
    double rad = (2*PI/(double)bundray_cnt)*(double)i;
    point p = {1000.0*sin(rad), 1000.0*cos(rad)};
    polygon.push_back(p);
  }
/*  for(auto out : polygon){
    cout << out.x <<"  "<< out.y <<"\n";
  };*/
  int v_end = v_baundray;
  int idx_p = 0;
  do{
    int h = outgoing_halfedges[v_baundray];
    while(opposite_halfedges[h] != -1){
      int h_opp = opposite_halfedges[h];
      if(h_opp%3 != 2) h = h_opp+1;
      else h = h_opp-2;
      //cout << h <<"\n";
    }
    int j = calc_v_tgt(F[h/3], h);
    bx(j) = polygon[idx_p].x;
    by(j) = polygon[idx_p].y;
    idx_p++;
    v_baundray = j;
    //cout <<"注目点"<< v_baundray <<"  "<< v_end <<"\n";
   // break;
  }while(v_baundray != v_end);
// 方程式を解く
  printf("方程式を解く\n");
  VectorXd x = laplacian.fullPivLu().solve(bx);
  VectorXd y = laplacian.fullPivLu().solve(by);

  ofstream file;
  file.open("project1.off", ios::out);
  file << "OFF\n";
  file << v_num <<"  "<< f_num <<"  0\n";
  for(int i=0; i<v_num; i++){
    file << x(i) <<"  "<< y(i) <<"  "<< z[i] <<"\n";
  }
  for(int i=0; i<f_num; i++){
    file << "3  " << F[i].v1 <<"  "<< F[i].v2 <<"  "<< F[i].v3 <<"\n";
  }
  file.close();
/*  SparseSolver solver;
  solver.compute(laplacian);

  Vector x(v_num);
  Vector y(v_num);
  x = solver.solve(bx);
  y = solver.solve(by);*/
  return 0;
}

