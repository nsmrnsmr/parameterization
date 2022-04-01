/*平面パラメータ化 Tuttle埋め込み
  実行までの例（Eigenのパスは自身のパスを指定してください）
  g++ param_spar.cpp -o param_spar -std=c++17 -I /usr/local/.../include/eigen3
  ./param_spar sample.off
*/
# include<iostream>
# include<cmath>
# include<cstdio>
# include<vector>
# include<algorithm>
# include<Eigen/Eigen>
# include<Eigen/Sparse>
# include<fstream>
# define PI 3.14159265358979323846
# define N 100
using SparseMatrixd = Eigen::SparseMatrix<double>;
using SparseSolver = Eigen::SparseLU<SparseMatrixd>;
using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using Tripletd = Eigen::Triplet<double>;
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
  for(int i=0; i<v_num; i++){
    char c[N] = {0};
    fgets(c, N, fl);
    double x, y, zz;
    sscanf(c, "%lf %lf %lf", &x, &y, &zz);
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
//メッシュのパラメータ化  //////////////////////////////////////////////////////////
  printf("メッシュのパラメータ化\n");
/*  MatrixXd laplacian(v_num, v_num); laplacian = MatrixXd::Zero(v_num, v_num);
  VectorXd bx(v_num); bx = VectorXd::Zero(v_num);
  VectorXd by(v_num); by = VectorXd::Zero(v_num);
*/
  SparseMatrixd laplacian(v_num, v_num);
  Vector bx = Vector::Zero(v_num);
  Vector by = Vector::Zero(v_num);

  vector<Tripletd> lap;
  int v_baundray = 0;
  int bundray_cnt = 0;
  for(int i=0; i<v_num; i++){ //点iに接続している頂点数の算出し、laplacianを求める
    int h = outgoing_halfedges[i]; //点iを支点とした半辺の添字
    int h_end = h; //終了条件（半時計周りに半辺を通して頂点を調べるためスタート地点を保存する）
    int cnt = 0; //繋がっている頂点数
    vector<Tripletd> lapval;
    do{
      int j = calc_v_tgt(F[h/3], h);
      lapval.push_back(Tripletd(i, j, -1)); // laplacian(i, j) = -1;
      cnt++;
      int h_prev = 0;
      if(h%3 != 0) h_prev = h-1;
      else h_prev = h+2;
      h = opposite_halfedges[h_prev];
      if(h == -1){ //境界辺に当たった場合の処理
        cnt = -1;
        bundray_cnt++;
        v_baundray = i;
        lap.push_back(Tripletd(i, i, 1)); //laplacian(i, i) = 1;
        break;
      }
    }while(h != h_end);
    if(cnt != -1){
      lapval.push_back(Tripletd(i, i, (double)cnt)); //laplacian(i, i) = cnt;
      lap.insert(lap.begin(), lapval.begin(), lapval.end());
    }
  }
  laplacian.setFromTriplets(lap.begin(), lap.end());

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
  SparseSolver solver;
  solver.compute(laplacian);
  VectorXd x = solver.solve(bx);
  VectorXd y = solver.solve(by);

  ofstream file;
  file.open("param_spar.off", ios::out);
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

