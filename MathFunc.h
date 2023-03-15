#pragma once

#include <math.h>
#include <vector>

// typedef double mf_matrix[4][4];
// typedef double mf_vector[4];

// namespace dgn
// {
// 	namespace math
// 	{
// 		class CVector3
// 		{
// 		public:
// 			enum enAxis { enAxisX = 0, enAxisY, enAxisZ, enAxisNum };
// 			// coordinates
// 			double x, y, z;

// 		public:
// 			// construction
// 			CVector3();
// 			CVector3(const double& x, const double& y, const double& z);
// 			CVector3(const double coordinate[enAxisNum]);
// 			CVector3(const CVector3& rVector);

// 		public:
// 			// operator overloading
// 			const double& operator[] (int i) const;
// 			double& operator[] (int i);

// 			CVector3& operator= (const CVector3& rVector);
// 			bool operator== (const CVector3& rVector) const;
// 			bool operator!= (const CVector3& rVector) const;

// 			CVector3 operator+ (const CVector3& v) const;
// 			CVector3 operator- (const CVector3& v) const;
// 			CVector3 operator* (double s) const;
// 			CVector3 operator/ (double s) const;
// 			CVector3 operator* (const CVector3& v) const;
// 			CVector3 operator/ (const CVector3& v) const;
// 			CVector3 operator- () const;

// 			CVector3& operator+= (const CVector3& v);
// 			CVector3& operator-= (const CVector3& v);
// 			CVector3& operator*= (double s);
// 			CVector3& operator/= (double s);
// 			CVector3& operator*= (const CVector3& v);
// 			CVector3& operator/= (const CVector3& v);

// 			// implement
// 			void GetCoord(double& x, double& y, double& z);
// 			void GetCoord(double coord[3]);
// 			bool isZero() const; // Returns true if this vector has length == 0.0
// 			bool isUnit() const; // Returns true if this vector has length == 1.0
// 			double length() const;
// 			double squaredLength() const;
// 			double dot(const CVector3& rVector) const;
// 			CVector3 cross(const CVector3& rkVector) const;
// 			double normalize();
// 			CVector3 normalizecross(const CVector3& rkVector) const;

// 			CVector3 minVec(const CVector3 &v) const;
// 			CVector3 maxVec(const CVector3 &v) const;

// 			// Linear interpolation
// 			inline CVector3 lerp(const CVector3& v, double alpha) const {
// 				return ( *this ) + ( v - *this ) * alpha;
// 			}

// 			inline double sum() const {
// 				return x + y + z;
// 			}

// 			inline double average() const {
// 				return sum() / 3.0;
// 			}

// 			inline static const CVector3& zero() { static CVector3 v(0.0, 0.0, 0.0); return v; }
// 			inline static const CVector3& unitX() { static CVector3 v(1.0, 0.0, 0.0); return v; }
// 			inline static const CVector3& unitY() { static CVector3 v(0.0, 1.0, 0.0); return v; }
// 			inline static const CVector3& unitZ() { static CVector3 v(0.0, 0.0, 1.0); return v; }
// 		};
// 	}
// }

class CDgnMathFunc
{
public:
	CDgnMathFunc() {};
	virtual ~CDgnMathFunc() {};

private:
	static double m_pi;
	static double m_trrad;
	static double m_trang;
	static double m_NormalZero;
	static double m_VertTol;
	static double m_AngleZero;

public:
	// static double Get_pi();
	// static double Get_trrad();
	// static double Get_trang();
	// static double Get_NormalZero();

	// static double GetAngle(double x1, double y1, double x2, double y2, bool bPosAngle = true);
	// static double CalcArea(const double& dDia);

	// static int math_ccw(double ax, double ay, double bx, double by, double cx, double cy);
	// static int math_ccw_AutoTol(double ax, double ay, double bx, double by, double cx, double cy);

	// static void mathMxIdentity(mf_matrix m);
	// static void mathMxMult(mf_matrix m1, mf_matrix m2, mf_matrix mout);
	// static void mathMxCopy(mf_matrix min, mf_matrix mout);

	// static void mathVecMult(mf_vector v, mf_matrix m, mf_vector vout);
	// static void mathVecCross(mf_vector v1, mf_vector v2, mf_vector vout);
	// static double mathVecDot(mf_vector v1, mf_vector v2);

	// static void mathRotateX(double angle, double& rx, double& ry, double& rz);  //Global-X축에 대한 회전, angle:회전각[deg], rx,ry,rz:회전대상 좌표
	// static void mathRotateY(double angle, double& rx, double& ry, double& rz);  //Global-Y축에 대한 회전, angle:회전각[deg], rx,ry,rz:회전대상 좌표
	// static void mathRotateZ(double angle, double& rx, double& ry, double& rz);  //Global-Z축에 대한 회전, angle:회전각[deg], rx,ry,rz:회전대상 좌표
	// static void mathRotate(double angle, double ux, double uy, double uz, double& rx, double& ry, double& rz);  //임의 축에 대한 회전, angle:회전각[deg], ux,uy,uz:회전축 벡터, rx,ry,rz:회전대상 좌표
	// static void mathRotate(double angle, double px, double py, double pz, double ux, double uy, double uz,
	// 	double& rx, double& ry, double& rz); //임의 축에 대한 회전, angle:회전각[deg], px,py,pz:회전축상의 한점, ux,uy,uz:회전축 벡터, rx,ry,rz:회전대상 좌표
	// static double mathDot(double xyz1[3], double xyz2[3], bool bNormal=true);
	// static double mathDot(double u1x, double u1y, double u2x, double u2y, bool bNormal = true);
	// static void mathCross(double xyz1[3], double xyz2[3], double xyzout[3]);
	// static double mathCrossAngleNormalize(double Vector1[3], double Vector2[3]);
	// static double mathCrossAngleNormalizeWithSign(double Vector1[3], double Vector2[3]);
	// static double mathCrossAngle(double Vector1[3], double Vector2[3]);                                 //두벡터 사이의 각도를 계산
	// static double mathCrossAngle2D(double u1x, double u1y, double u2x, double u2y);                     //normalize  된 두 2차원벡터 사이의 각도를 계산
	// static double mathCrossAngle2D(double Vector1[3], double Vector2[3]);                                 //두 2차원벡터 사이의 각도를 계산 / +-180도 범위로.
	// static double mathLength(double dx, double dy, double dz = 0.);                                       //1점과 원점과의 거리(2,3차원)
	// static double mathLength(double dxi, double dyi, double dxj, double dyj);                           //2점간의 거리(2차원)
	// static double mathLength(double dxi, double dyi, double dzi, double dxj, double dyj, double dzj);   //2점간의 거리(3차원)
	// static double mathLength(double xyz1[3], double xyz2[3]);   //2점간의 거리(3차원)
	// static double mathArea(double xyz1[3], double xyz2[3], double xyzout[3]);                           //삼각형 면적
	// static double mathArea2d(double x1, double y1, double x2, double y2, double x3, double y3);         //삼각형 면적(2차원)
	// static double mathArea2d(int nPoint, double x[], double y[]);
	// static double mathVolume(int nVertex, double xyz[][3]);                                             // 4,6,8 개의 꼭지점을 가진 고체요소의 체적
	// static void mathNormal(double Vector1[3], double Vector2[3], double VectorN[3]);                    //원점기준의 2선의 법선
	// static void mathNormal(double Vector1[3], double Vector2[3], double Vector3[3], double VectorN[3]);//Vector1기준의 2선의 법선
	// static bool mathNormalize(double xyz[3], double xyzn[3]);                                           //단위벡터로 변환
	// static bool mathNormalize(double dx, double dy, double dz, double& dxn, double& dyn, double& dzn);  //단위벡터로 변환(3d)
	// static bool mathNormalize(double dx, double dy, double& dxn, double& dyn);                          //단위벡터로 변환(2d)
	// static void mathPLCrossPoint(double a, double b, double c, double l, double m, double n,
	// 	double x, double y, double z, double& rpx, double& rpy, double& rpz);  //선과 점의 수직교점, a,b,c:선상의 1점, l,m,n:선의 방향벡터, x,y,z:점의 좌표, rpx,rpy,rpz:교점의 좌표
	// static void mathPLCrossPoint(double a, double b, double c, double l, double m, double n,
	// 	double x, double y, double z, double& rpx, double& rpy, double& rpz, double& t);  //선과 점의 수직교점, a,b,c:선상의 1점, l,m,n:선의 방향벡터, x,y,z:점의 좌표, rpx,rpy,rpz:교점의 좌표
	// static void mathVec2Ang(double ux, double uy, double uz, double& rx, double& ry);  //벡터를 각도로 변환(Z축을 Y축으로 ry, X축으로 rx 회전시키면 ux,uy,uz위치가 됨)
	// static bool mathPlaneEquation(double p1[3], double p2[3], double p3[3], double& a, double& b, double& c, double & h);
	// static double mathDistanceToPoint(double p1[3], double p2[3]);
	// static double mathDistanceToLine(double line_i[3], double line_j[3], double point[3]); //점과 선분과의 거리
	// static double mathDistanceToLine2D(double line_i_x, double line_i_y, double line_j_x, double line_j_y, double point_x, double point_y);
	// static double mathDistanceToPlane(double a, double b, double c, double h, double x, double y, double z); // 점과 면의 수직거리
	// static double mathDistanceToPlaneWithSign(double a, double b, double c, double h, double x, double y, double z); // 점과 면의 수직거리 부호 포함

	// /////////////////////////////////////////////////
	// // Add by sshan(090708)
	// // 평면에 점의 포함여부 반환
	// // p1, p2, p3 : Plane을 형성할 포인터
	// // point : 포함여부 판단할 점
	// static bool mathIncludePointInPlane(double p1[3], double p2[3], double p3[3], double point[3], double Tol = 0.);

	// // Add by ZINU.('02.09.26). 점과 직선사이의 수직 거리
	// static double mathDistanceFromIntersectPointToLine(double line_i[3], double line_j[3], double point[3]);
	// static bool mathIntersectLine2D(double p1i_x, double p1i_y, double p1j_x, double p1j_y, double p2i_x, double p2i_y, double p2j_x, double p2j_y, double& pInts_x, double& pInts_y, bool bCheckOnLine = false); // 두개의 선분의 교점
	// static bool mathIntersectLinePMCv(double p1i_x, double p1i_y, double p1j_x, double p1j_y, double p2i_x, double p2i_y, double p2j_x, double p2j_y, double& pInts_x, double& pInts_y);
	// // 점과 직선의 교점
	// static bool mathIntersectPointToLine(double line_i[3], double line_j[3], double point[3], double rpInt[3]);
	// // 직선과 평면과의 교점계산
	// static bool mathIntersectPointToPlane(double LineVector[3], double LinePoint[3], double PlaneNormal[3], double PlanePoint[3], double pInts[3]);
	// // 직선과 평면과의 교점계산, dDistance:교점이 직선내에 있을경우에는 0. 아니면 교점과의 거리
	// static bool mathIntersectPointToPlane(double LinePoint1[3], double LinePoint2[3], double PlaneNormal[3], double PlanePoint[3], double pInts[3], double& dDistance);
	// static double mathDistanceLine(double pl1_i[3], double pl1_j[3], double pl2_i[3], double pl2_j[3]); // 두개의 선분의 최단 거리
	// static bool mathIntersectLine(double pl1_i[3], double pl1_j[3], double pl2_i[3], double pl2_j[3], double Tolerance, double& dblDistance, double pInts[3]); // 두개의 선분의 교점
	// static bool mathIntersectLine2(double pl1_i[3], double pl1_j[3], double pl2_i[3], double pl2_j[3], double Tolerance, double& dblDistance, double pInts[3]); // 두개의 선분 내의 교점
	// // 임의 좌표계로 좌표변환(ux,uy,uz는 단위 벡터, 단 원점이동은 없음)
	// static void mathTranUCS(const double ux[3], const double uy[3], const double uz[3], const int nData, double coor[][3]);
	// static void mathGCS2UCS(double& x, double& y, double& z, double ucs[3][3]);
	// static void mathUCS2GCS(double& x, double& y, double& z, double ucs[3][3]);
	// static bool mathPolyCentroid(int n, double x[], double y[], double& xCentroid, double& yCentroid, double& area);
	// // 직교 좌표계의 한 점을 원통 좌표계로 바꾸는 함수(x, y, z) => (radius, angle, height)
	// static void mathConvertToCylinderCoordi(double node_xyz[3], double cyn_xyz[3], double org_xyz[3], double rot_xyz[3], double pol_xyz[3]);
	// // 대상 좌표    // 변환된 원통 좌표    // 극점좌표    // 회전축 위의 좌표// 극축 위의 좌표
	// static bool project_on_cylinder_normal(const double node_xyz[3], const double radius, const double height, const double tolerance, double projected_node_xyz[3]);
	// static bool project_on_cylinder_vector(const double node_xyz[3], const double vector_xyz[3], const double radius, const double height, const double tolerance, double projected_node_xyz[3]);
	// static bool project_on_cone_normal(const double node_xyz[3], const double radius1, const double radius2, const double height, const double tolerance, double projected_node_xyz[3]);
	// static bool project_on_cone_vector(const double node_xyz[3], const double vector_xyz[3], const double radius1, const double radius2, const double height, const double tolerance, double projected_node_xyz[3]);
	// static bool project_on_sphere_normal(const double node_xyz[3], const double radius, const double tolerance, double projected_node_xyz[3]);
	// static bool project_on_sphere_vector(const double node_xyz[3], const double vector_xyz[3], const double radius, const double tolerance, double projected_node_xyz[3]);
	// static bool project_on_ellipsoid_normal(const double node_xyz[3], const double radius, const double height, const double tolerance, double projected_node_xyz[3]);
	// static bool project_on_ellipsoid_vector(const double node_xyz[3], const double vector_xyz[3], const double radius, const double height, const double tolerance, double projected_node_xyz[3]);
	// static bool project_on_plane_normal(const double node_xyz[3], const double plane_xyz1[3], const double plane_xyz2[3], const double plane_xyz3[3], const double tolerance, double projected_node_xyz[3]);
	// static bool project_on_plane_vector(const double node_xyz[3], const double vector_xyz[3], const double plane_xyz1[3], const double plane_xyz2[3], const double plane_xyz3[3], double projected_node_xyz[3], double tolerance);
	// static bool project_on_plane_vector(const double node_xyz[3], const double vector_xyz[3], const double plane_xyz1[3], const double plane_xyz2[3], const double plane_xyz3[3], double projected_node_xyz[3]);
	// static bool project_on_plane_vector_for_CuttingPlane(const double node_xyz[3], const double vector_xyz[3], const double plane_xyz1[3], const double plane_xyz2[3], const double plane_xyz3[3], double projected_node_xyz[3]);
	// static bool project_on_quad_curve_normal(const double curve_xyz_1[3], const double curve_xyz_2[3], const double curve_xyz_3[3], const double point[3], const double tolerance, double project_point[3]);
	// static bool project_on_quad_curve_vector(const double curve_xyz_1[3], const double curve_xyz_2[3], const double curve_xyz_3[3], const double point[3], const double direction_vector[3], const double tolerance, double project_point[3]);
	// static bool check_quad_curve(const double curve_xyz_1[3], const double curve_xyz_2[3], const double curve_xyz_3[3], const double ZERO);

	// static double mathInterpolate(double dLeft, double dRight, double dDistRatioFromLeft);
	// static double mathInterpolate(double dLeftBot, double dRightBot, double dLeftTop, double dRightTop, double dDistRatioFromLeft, double dDistRatioFromBot);
	// static double InterpolateLinear(const double& dX1, const double& dX2, const double& dY1, const double& dY2, const double& dX);
	// static double InterpolateLinear(const std::vector<std::pair<double, double>>& vdXY, const double& dX);

	// // SHIN (2006.2.3) 2차원 도형계산을 위해 추가함

	// static void mathCrossPoint2D(double line_i[2], double line_j[2], double point[2], double point_cross[2]);
	// // 2차원상의 점과 직선사이의 수직 거리(nState : 해당 점과의 직교점위 위치 (0:i밖에 있을때 1:선분내에 있을때 2:j밖에 있을때))
	// static double mathDistanceFromIntersectPointToLine2D(double line_i[2], double line_j[2], double point[2], int& nState);
	// //normalize  된 두 2차원벡터 사이의 각도를 계산 +-180도 범위의 값을 가짐
	// static double mathCrossAngle2DSign(double u1x, double u1y, double u2x, double u2y);
	// // 두 2차원 vector의 외적을 구하여 넘겨줌
	// static double mathCross2D(double vector1[2], double vector2[2]);
	// // 2차원 vector의 크기를 단위vector로 변환
	// //   vector  : 변환할 벡터
	// //   vectorn : Normalize된 벡터를 넘겨받을 변수
	// static bool mathNormalize2D(double vector[2], double vectorn[2]);
	// // 2차원 좌표 p1과 p2의 값을 서로 바꾸어줌
	// static void mathSwap2D(double p1[2], double p2[2]);
	// // 2차원 Line내에 해당좌표가 포함되는지 여부를 판단함
	// //   bound1  : Line의 구성 좌표
	// //   bound2  : Line의 구성 좌표
	// //   targetPt: 포함여부를 확인할 좌표
	// //   isOnLine: bound좌표와 targetPt좌표가 일치할경우 포함으로 볼지 여부
	// //   return  : true=포함   false=미포함
	// static bool mathIsPointOfLine2D(double bound1[2], double bound2[2], double targetPt[2], bool isOnLine);
	// // 해당 좌표(p1)가 2차원 polyLine내에 들어가는지 여부를 판단함
	// //   p1      : 포함여부를 확인할 좌표
	// //   nData   : polyLine의 좌표 갯수
	// //   polyLine: 비교할 2차원 좌표의 배열
	// //   bIncludeOutLine : 해당좌표가 Line 및 꼭지점에 일치할경우 포함으로 볼지 여부
	// //   return  : true=포함   false=미포함
	// // ※참고문헌(기하 알고리즘 P1124)
	// static bool mathIsInsidePoint2D(double p1[2], const int nData, double polyLine[][2], bool bIncludeOutLine);
	// // 두 선분이 교차하는지 여부를 판단함
	// //   p1Org : 첫번째 선분의 구성좌표
	// //   p2Org : 첫번째 선분의 구성좌표
	// //   p3Org : 두번째 선분의 구성좌표
	// //   p4Org : 두번째 선분의 구성좌표
	// //   return: 1=교차하는 경우  -1=교차하지 않는경우  0=한선분의 끝점이 다른선분에 포함되는 경우
	// // ※참고문헌(기하 알고리즘 P1118)
	// static int  mathIntersect_ccw2D(double p1Org[2], double p2Org[2], double p3Org[2], double p4Org[2]);
	// static int  mathIntersect_ccw2D_Tol(double p1Org[2], double p2Org[2], double p3Org[2], double p4Org[2]);

	// // 직선과 직선의 교점 구하기
	// //   line1 : Line1을 구성하는 2개의 좌표점 (line1[1][0] 2번째 좌표의 x값, line1[1][1] 2번째 좌표의 y값)
	// //   line2 : Line2을 구성하는 2개의 좌표점
	// //   cross : 교차점의 좌표
	// //   return: 0:평행  1:교차
	// static int mathLineLineCross2D(double line1[2][2], double line2[2][2], double cross[2]);
	// // 직선과 선분의 교점 구하기
	// //   line1 : Line을 구성하는 2개의 좌표점 (line1[1][0] 2번째 좌표의 x값, line1[1][1] 2번째 좌표의 y값)
	// //   sttP  : 선분을 구성하는 좌표점
	// //   endP  : 선분을 구성하는 좌표점
	// //   cross : 교차점의 좌표
	// //   return: 0:교차점이 선분내에 없을때  1:교차점이 선분내에 있을때
	// static int mathLineSegCross2D(double line[2][2], double sttP[2], double endP[2], double cross[2]);

	// // 2차원 polyline을 dOffset거리만큼 도형의 바깥쪽, 안쪽으로 크기를 변경함
	// //   dOffset : Offset거리(-:도형의 안쪽으로 Offset,  +:도형의 바깥쪽으로 Offset)
	// //   nData   : polyLine의 좌표 갯수
	// // < polyLine: 입력된 polyLine을 Offset시켜서 넘겨줌
	// //   nRotType: polyLine의 회전방향(0:Auto,  1:반시계방향,  2:시계방향)
	// static bool mathOffsetOfPolyline(double dOffset, const int nData, double polyLine[][2], int nRotType = 0);

	// // 3점으로 부터 원의 정보를 계산
	// //   P1, P2, P3 : 호를 구성하는 3점(1,2,3순서)
	// // < CenterP    : 원의 중심점
	// // < Radius     : 반지름
	// //   RETURN     : -1:3점이 모두 같은 경우, 0:3점이 직선인경우, 1:계산가능한 경우
	// static int  mathCircleForm3Point(double p1[2], double p2[2], double p3[2], double CenterP[2], double& dRadius);
	// // 3점으로 부터 호의 정보를 계산
	// //   P1, P2, P3 : 호를 구성하는 3점(1,2,3순서)
	// // < CenterP    : 원의 중심점
	// // < Radius     : 반지름
	// // < StartAngle : 시작각도
	// // < SweepAngle : 내부각도(반시계방향 +)
	// //   RETURN     : -1:3점이 모두 같은 경우, 0:3점이 직선인경우, 1:계산가능한 경우
	// static int  mathArcForm3Point(double p1[2], double p2[2], double p3[2], double CenterP[2], double& dRadius, double& dStartAngle, double& dSweepAngle);
	// static int  mathArcForm3Point(double p1[2], double p2[2], double p3[2], double CenterP[2], double& dRadius, double& dStartAngle, double& dSweepAngle, double& dBulge);

	// static bool   GetCenter_Bulge(double line[2][2], double bulge, double centerP[2]);
	// static double GetBulge(double line[2][2], double radius, bool inside);
	// static int    GetRotationType(const int nData, double polyLine[][2]);//polyLine의 회전방향을 넘겨줌(1:반시계방향,  2:시계방향)

	// // 해당평면을 기준으로 만들어진 임의 형상의 기둥(상면은 평면이어야함)의 부피 및 중심점을 계산하여 넘겨줌
	// // < pointCenV, dVolm : 계산된 중심점과 부피
	// //   nData, PontUnitList : 기둥 상면의 3D좌표 및 점의 개수
	// //   nAxis : 계산할 평면의 방향 (1:xy평면, 2:yz평면, 3:xz평면)
	// static bool Get_VolmCenter(double pointCenV[3], double& dVolm, const int nData, double PontUnitList[][3], int nAxis = 1);

	// // 두 Polyline의 영역 중첩여부 체크(교차 되거나 내부포함 겹칠경우 true, 겹치는 곳이 없을 경우 false) - Test전 사용보류(by SHIN 2010.05.31)
	// static bool IsOverLap_Polyline(const int nData1, double polyLine1[][2], const int nData2, double polyLine2[][2]);

	// // 유한평면과 직선과의 교점
	// static bool IntersectTriangle(double RayOrigin[3], double RayDirection[3], double V0[3], double V1[3], double V2[3], double cross[3]);
	// static bool IntersectTriangle(const dgn::math::CVector3& RayOrigin, const dgn::math::CVector3& RayDirection,
	// 	const dgn::math::CVector3& V0, const dgn::math::CVector3& V1, const dgn::math::CVector3& V2, dgn::math::CVector3& cross);

	static double mathSin(double x);
	static double mathCos(double x);
	static double mathTan(double x);
	static double mathCot(double x);
	static double mathAsin(double x);
	static double mathAcos(double x);
	static double mathSqrt(double x);
	static int    mathAbs(int x);
	static double mathAbs(double x);

	static double mathRoundOff(double dVal, int N);
	static double mathRoundUp(double dVal, int N);
	static double mathRoundDown(double dVal, int N);

	static bool   mathFrameLocalVector(double Coor_i[3], double Cood_j[3], double Angle, double FrameLocalVector[3][3]);
	static bool   mathPlaneLocalVector(int NumNode, double Coor[][3], double dAngle, double PlaneLocalVector[3][3]);

protected:
	static double Solid_4(double xyz[][3]);
	static double Solid_6(double xyz[][3]);
	static double Solid_8(double xyz[][3]);

public:
	// static int RealToInt(const double dVal, double dTol = 1e-6);

// public:
// 	static double mathMomentOfInertiaRectangle(double dB, double dH, double de);
// 	static double mathMomentOfInertiaTriangle(double dB, double dH, double de);
// 	static double mathMomentOfInertiaCircle(double dr);

private:
	static bool IsNormalZero(const double& dVal);
};
