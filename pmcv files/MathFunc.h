#pragma once

#include <math.h>
#include "HeaderPre.h"
#include <vector>

typedef double mf_matrix[4][4];
typedef double mf_vector[4];

namespace dgn
{
	namespace math
	{
		class CVector3
		{
		public:
			enum enAxis { enAxisX = 0, enAxisY, enAxisZ, enAxisNum };
			// coordinates
			double x, y, z;

		public:
			// construction
			CVector3();
			CVector3(const double& x, const double& y, const double& z);
			CVector3(const double coordinate[enAxisNum]);
			CVector3(const CVector3& rVector);

		public:
			// operator overloading
			const double& operator[] (int i) const;
			double& operator[] (int i);

			CVector3& operator= (const CVector3& rVector);
			bool operator== (const CVector3& rVector) const;
			bool operator!= (const CVector3& rVector) const;

			CVector3 operator+ (const CVector3& v) const;
			CVector3 operator- (const CVector3& v) const;
			CVector3 operator* (double s) const;
			CVector3 operator/ (double s) const;
			CVector3 operator* (const CVector3& v) const;
			CVector3 operator/ (const CVector3& v) const;
			CVector3 operator- () const;

			CVector3& operator+= (const CVector3& v);
			CVector3& operator-= (const CVector3& v);
			CVector3& operator*= (double s);
			CVector3& operator/= (double s);
			CVector3& operator*= (const CVector3& v);
			CVector3& operator/= (const CVector3& v);

			// implement
			void GetCoord(double& x, double& y, double& z);
			void GetCoord(double coord[3]);
			bool isZero() const; // Returns true if this vector has length == 0.0
			bool isUnit() const; // Returns true if this vector has length == 1.0
			double length() const;
			double squaredLength() const;
			double dot(const CVector3& rVector) const;
			CVector3 cross(const CVector3& rkVector) const;
			double normalize();
			CVector3 normalizecross(const CVector3& rkVector) const;

			CVector3 minVec(const CVector3 &v) const;
			CVector3 maxVec(const CVector3 &v) const;

			// Linear interpolation
			inline CVector3 lerp(const CVector3& v, double alpha) const {
				return ( *this ) + ( v - *this ) * alpha;
			}

			inline double sum() const {
				return x + y + z;
			}

			inline double average() const {
				return sum() / 3.0;
			}

			inline static const CVector3& zero() { static CVector3 v(0.0, 0.0, 0.0); return v; }
			inline static const CVector3& unitX() { static CVector3 v(1.0, 0.0, 0.0); return v; }
			inline static const CVector3& unitY() { static CVector3 v(0.0, 1.0, 0.0); return v; }
			inline static const CVector3& unitZ() { static CVector3 v(0.0, 0.0, 1.0); return v; }
		};
	}
}

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
	static double Get_pi();
	static double Get_trrad();
	static double Get_trang();
	static double Get_NormalZero();

	static double GetAngle(double x1, double y1, double x2, double y2, bool bPosAngle = true);
	static double CalcArea(const double& dDia);

	static int math_ccw(double ax, double ay, double bx, double by, double cx, double cy);
	static int math_ccw_AutoTol(double ax, double ay, double bx, double by, double cx, double cy);

	static void mathMxIdentity(mf_matrix m);
	static void mathMxMult(mf_matrix m1, mf_matrix m2, mf_matrix mout);
	static void mathMxCopy(mf_matrix min, mf_matrix mout);

	static void mathVecMult(mf_vector v, mf_matrix m, mf_vector vout);
	static void mathVecCross(mf_vector v1, mf_vector v2, mf_vector vout);
	static double mathVecDot(mf_vector v1, mf_vector v2);

	static void mathRotateX(double angle, double& rx, double& ry, double& rz);  //Global-X�࿡ ���� ȸ��, angle:ȸ����[deg], rx,ry,rz:ȸ����� ��ǥ
	static void mathRotateY(double angle, double& rx, double& ry, double& rz);  //Global-Y�࿡ ���� ȸ��, angle:ȸ����[deg], rx,ry,rz:ȸ����� ��ǥ
	static void mathRotateZ(double angle, double& rx, double& ry, double& rz);  //Global-Z�࿡ ���� ȸ��, angle:ȸ����[deg], rx,ry,rz:ȸ����� ��ǥ
	static void mathRotate(double angle, double ux, double uy, double uz, double& rx, double& ry, double& rz);  //���� �࿡ ���� ȸ��, angle:ȸ����[deg], ux,uy,uz:ȸ���� ����, rx,ry,rz:ȸ����� ��ǥ
	static void mathRotate(double angle, double px, double py, double pz, double ux, double uy, double uz,
		double& rx, double& ry, double& rz); //���� �࿡ ���� ȸ��, angle:ȸ����[deg], px,py,pz:ȸ������� ����, ux,uy,uz:ȸ���� ����, rx,ry,rz:ȸ����� ��ǥ
	static double mathDot(double xyz1[3], double xyz2[3], bool bNormal=true);
	static double mathDot(double u1x, double u1y, double u2x, double u2y, bool bNormal = true);
	static void mathCross(double xyz1[3], double xyz2[3], double xyzout[3]);
	static double mathCrossAngleNormalize(double Vector1[3], double Vector2[3]);
	static double mathCrossAngleNormalizeWithSign(double Vector1[3], double Vector2[3]);
	static double mathCrossAngle(double Vector1[3], double Vector2[3]);                                 //�κ��� ������ ������ ���
	static double mathCrossAngle2D(double u1x, double u1y, double u2x, double u2y);                     //normalize  �� �� 2�������� ������ ������ ���
	static double mathCrossAngle2D(double Vector1[3], double Vector2[3]);                                 //�� 2�������� ������ ������ ��� / +-180�� ������.
	static double mathLength(double dx, double dy, double dz = 0.);                                       //1���� �������� �Ÿ�(2,3����)
	static double mathLength(double dxi, double dyi, double dxj, double dyj);                           //2������ �Ÿ�(2����)
	static double mathLength(double dxi, double dyi, double dzi, double dxj, double dyj, double dzj);   //2������ �Ÿ�(3����)
	static double mathLength(double xyz1[3], double xyz2[3]);   //2������ �Ÿ�(3����)
	static double mathArea(double xyz1[3], double xyz2[3], double xyzout[3]);                           //�ﰢ�� ����
	static double mathArea2d(double x1, double y1, double x2, double y2, double x3, double y3);         //�ﰢ�� ����(2����)
	static double mathArea2d(int nPoint, double x[], double y[]);
	static double mathVolume(int nVertex, double xyz[][3]);                                             // 4,6,8 ���� �������� ���� ��ü����� ü��
	static void mathNormal(double Vector1[3], double Vector2[3], double VectorN[3]);                    //���������� 2���� ����
	static void mathNormal(double Vector1[3], double Vector2[3], double Vector3[3], double VectorN[3]);//Vector1������ 2���� ����
	static bool mathNormalize(double xyz[3], double xyzn[3]);                                           //�������ͷ� ��ȯ
	static bool mathNormalize(double dx, double dy, double dz, double& dxn, double& dyn, double& dzn);  //�������ͷ� ��ȯ(3d)
	static bool mathNormalize(double dx, double dy, double& dxn, double& dyn);                          //�������ͷ� ��ȯ(2d)
	static void mathPLCrossPoint(double a, double b, double c, double l, double m, double n,
		double x, double y, double z, double& rpx, double& rpy, double& rpz);  //���� ���� ��������, a,b,c:������ 1��, l,m,n:���� ���⺤��, x,y,z:���� ��ǥ, rpx,rpy,rpz:������ ��ǥ
	static void mathPLCrossPoint(double a, double b, double c, double l, double m, double n,
		double x, double y, double z, double& rpx, double& rpy, double& rpz, double& t);  //���� ���� ��������, a,b,c:������ 1��, l,m,n:���� ���⺤��, x,y,z:���� ��ǥ, rpx,rpy,rpz:������ ��ǥ
	static void mathVec2Ang(double ux, double uy, double uz, double& rx, double& ry);  //���͸� ������ ��ȯ(Z���� Y������ ry, X������ rx ȸ����Ű�� ux,uy,uz��ġ�� ��)
	static bool mathPlaneEquation(double p1[3], double p2[3], double p3[3], double& a, double& b, double& c, double & h);
	static double mathDistanceToPoint(double p1[3], double p2[3]);
	static double mathDistanceToLine(double line_i[3], double line_j[3], double point[3]); //���� ���а��� �Ÿ�
	static double mathDistanceToLine2D(double line_i_x, double line_i_y, double line_j_x, double line_j_y, double point_x, double point_y);
	static double mathDistanceToPlane(double a, double b, double c, double h, double x, double y, double z); // ���� ���� �����Ÿ�
	static double mathDistanceToPlaneWithSign(double a, double b, double c, double h, double x, double y, double z); // ���� ���� �����Ÿ� ��ȣ ����

	/////////////////////////////////////////////////
	// Add by sshan(090708)
	// ��鿡 ���� ���Կ��� ��ȯ
	// p1, p2, p3 : Plane�� ������ ������
	// point : ���Կ��� �Ǵ��� ��
	static bool mathIncludePointInPlane(double p1[3], double p2[3], double p3[3], double point[3], double Tol = 0.);

	// Add by ZINU.('02.09.26). ���� ���������� ���� �Ÿ�
	static double mathDistanceFromIntersectPointToLine(double line_i[3], double line_j[3], double point[3]);
	static bool mathIntersectLine2D(double p1i_x, double p1i_y, double p1j_x, double p1j_y, double p2i_x, double p2i_y, double p2j_x, double p2j_y, double& pInts_x, double& pInts_y, bool bCheckOnLine = false); // �ΰ��� ������ ����
	static bool mathIntersectLinePMCv(double p1i_x, double p1i_y, double p1j_x, double p1j_y, double p2i_x, double p2i_y, double p2j_x, double p2j_y, double& pInts_x, double& pInts_y);
	// ���� ������ ����
	static bool mathIntersectPointToLine(double line_i[3], double line_j[3], double point[3], double rpInt[3]);
	// ������ ������ �������
	static bool mathIntersectPointToPlane(double LineVector[3], double LinePoint[3], double PlaneNormal[3], double PlanePoint[3], double pInts[3]);
	// ������ ������ �������, dDistance:������ �������� ������쿡�� 0. �ƴϸ� �������� �Ÿ�
	static bool mathIntersectPointToPlane(double LinePoint1[3], double LinePoint2[3], double PlaneNormal[3], double PlanePoint[3], double pInts[3], double& dDistance);
	static double mathDistanceLine(double pl1_i[3], double pl1_j[3], double pl2_i[3], double pl2_j[3]); // �ΰ��� ������ �ִ� �Ÿ�
	static bool mathIntersectLine(double pl1_i[3], double pl1_j[3], double pl2_i[3], double pl2_j[3], double Tolerance, double& dblDistance, double pInts[3]); // �ΰ��� ������ ����
	static bool mathIntersectLine2(double pl1_i[3], double pl1_j[3], double pl2_i[3], double pl2_j[3], double Tolerance, double& dblDistance, double pInts[3]); // �ΰ��� ���� ���� ����
	// ���� ��ǥ��� ��ǥ��ȯ(ux,uy,uz�� ���� ����, �� �����̵��� ����)
	static void mathTranUCS(const double ux[3], const double uy[3], const double uz[3], const int nData, double coor[][3]);
	static void mathGCS2UCS(double& x, double& y, double& z, double ucs[3][3]);
	static void mathUCS2GCS(double& x, double& y, double& z, double ucs[3][3]);
	static bool mathPolyCentroid(int n, double x[], double y[], double& xCentroid, double& yCentroid, double& area);
	// ���� ��ǥ���� �� ���� ���� ��ǥ��� �ٲٴ� �Լ�(x, y, z) => (radius, angle, height)
	static void mathConvertToCylinderCoordi(double node_xyz[3], double cyn_xyz[3], double org_xyz[3], double rot_xyz[3], double pol_xyz[3]);
	// ��� ��ǥ    // ��ȯ�� ���� ��ǥ    // ������ǥ    // ȸ���� ���� ��ǥ// ���� ���� ��ǥ
	static bool project_on_cylinder_normal(const double node_xyz[3], const double radius, const double height, const double tolerance, double projected_node_xyz[3]);
	static bool project_on_cylinder_vector(const double node_xyz[3], const double vector_xyz[3], const double radius, const double height, const double tolerance, double projected_node_xyz[3]);
	static bool project_on_cone_normal(const double node_xyz[3], const double radius1, const double radius2, const double height, const double tolerance, double projected_node_xyz[3]);
	static bool project_on_cone_vector(const double node_xyz[3], const double vector_xyz[3], const double radius1, const double radius2, const double height, const double tolerance, double projected_node_xyz[3]);
	static bool project_on_sphere_normal(const double node_xyz[3], const double radius, const double tolerance, double projected_node_xyz[3]);
	static bool project_on_sphere_vector(const double node_xyz[3], const double vector_xyz[3], const double radius, const double tolerance, double projected_node_xyz[3]);
	static bool project_on_ellipsoid_normal(const double node_xyz[3], const double radius, const double height, const double tolerance, double projected_node_xyz[3]);
	static bool project_on_ellipsoid_vector(const double node_xyz[3], const double vector_xyz[3], const double radius, const double height, const double tolerance, double projected_node_xyz[3]);
	static bool project_on_plane_normal(const double node_xyz[3], const double plane_xyz1[3], const double plane_xyz2[3], const double plane_xyz3[3], const double tolerance, double projected_node_xyz[3]);
	static bool project_on_plane_vector(const double node_xyz[3], const double vector_xyz[3], const double plane_xyz1[3], const double plane_xyz2[3], const double plane_xyz3[3], double projected_node_xyz[3], double tolerance);
	static bool project_on_plane_vector(const double node_xyz[3], const double vector_xyz[3], const double plane_xyz1[3], const double plane_xyz2[3], const double plane_xyz3[3], double projected_node_xyz[3]);
	static bool project_on_plane_vector_for_CuttingPlane(const double node_xyz[3], const double vector_xyz[3], const double plane_xyz1[3], const double plane_xyz2[3], const double plane_xyz3[3], double projected_node_xyz[3]);
	static bool project_on_quad_curve_normal(const double curve_xyz_1[3], const double curve_xyz_2[3], const double curve_xyz_3[3], const double point[3], const double tolerance, double project_point[3]);
	static bool project_on_quad_curve_vector(const double curve_xyz_1[3], const double curve_xyz_2[3], const double curve_xyz_3[3], const double point[3], const double direction_vector[3], const double tolerance, double project_point[3]);
	static bool check_quad_curve(const double curve_xyz_1[3], const double curve_xyz_2[3], const double curve_xyz_3[3], const double ZERO);

	static double mathInterpolate(double dLeft, double dRight, double dDistRatioFromLeft);
	static double mathInterpolate(double dLeftBot, double dRightBot, double dLeftTop, double dRightTop, double dDistRatioFromLeft, double dDistRatioFromBot);
	static double InterpolateLinear(const double& dX1, const double& dX2, const double& dY1, const double& dY2, const double& dX);
	static double InterpolateLinear(const std::vector<std::pair<double, double>>& vdXY, const double& dX);

	// SHIN (2006.2.3) 2���� ��������� ���� �߰���

	static void mathCrossPoint2D(double line_i[2], double line_j[2], double point[2], double point_cross[2]);
	// 2�������� ���� ���������� ���� �Ÿ�(nState : �ش� ������ �������� ��ġ (0:i�ۿ� ������ 1:���г��� ������ 2:j�ۿ� ������))
	static double mathDistanceFromIntersectPointToLine2D(double line_i[2], double line_j[2], double point[2], int& nState);
	//normalize  �� �� 2�������� ������ ������ ��� +-180�� ������ ���� ����
	static double mathCrossAngle2DSign(double u1x, double u1y, double u2x, double u2y);
	// �� 2���� vector�� ������ ���Ͽ� �Ѱ���
	static double mathCross2D(double vector1[2], double vector2[2]);
	// 2���� vector�� ũ�⸦ ����vector�� ��ȯ
	//   vector  : ��ȯ�� ����
	//   vectorn : Normalize�� ���͸� �Ѱܹ��� ����
	static bool mathNormalize2D(double vector[2], double vectorn[2]);
	// 2���� ��ǥ p1�� p2�� ���� ���� �ٲپ���
	static void mathSwap2D(double p1[2], double p2[2]);
	// 2���� Line���� �ش���ǥ�� ���ԵǴ��� ���θ� �Ǵ���
	//   bound1  : Line�� ���� ��ǥ
	//   bound2  : Line�� ���� ��ǥ
	//   targetPt: ���Կ��θ� Ȯ���� ��ǥ
	//   isOnLine: bound��ǥ�� targetPt��ǥ�� ��ġ�Ұ�� �������� ���� ����
	//   return  : true=����   false=������
	static bool mathIsPointOfLine2D(double bound1[2], double bound2[2], double targetPt[2], bool isOnLine);
	// �ش� ��ǥ(p1)�� 2���� polyLine���� ������ ���θ� �Ǵ���
	//   p1      : ���Կ��θ� Ȯ���� ��ǥ
	//   nData   : polyLine�� ��ǥ ����
	//   polyLine: ���� 2���� ��ǥ�� �迭
	//   bIncludeOutLine : �ش���ǥ�� Line �� �������� ��ġ�Ұ�� �������� ���� ����
	//   return  : true=����   false=������
	// ��������(���� �˰��� P1124)
	static bool mathIsInsidePoint2D(double p1[2], const int nData, double polyLine[][2], bool bIncludeOutLine);
	// �� ������ �����ϴ��� ���θ� �Ǵ���
	//   p1Org : ù��° ������ ������ǥ
	//   p2Org : ù��° ������ ������ǥ
	//   p3Org : �ι�° ������ ������ǥ
	//   p4Org : �ι�° ������ ������ǥ
	//   return: 1=�����ϴ� ���  -1=�������� �ʴ°��  0=�Ѽ����� ������ �ٸ����п� ���ԵǴ� ���
	// ��������(���� �˰��� P1118)
	static int  mathIntersect_ccw2D(double p1Org[2], double p2Org[2], double p3Org[2], double p4Org[2]);
	static int  mathIntersect_ccw2D_Tol(double p1Org[2], double p2Org[2], double p3Org[2], double p4Org[2]);

	// ������ ������ ���� ���ϱ�
	//   line1 : Line1�� �����ϴ� 2���� ��ǥ�� (line1[1][0] 2��° ��ǥ�� x��, line1[1][1] 2��° ��ǥ�� y��)
	//   line2 : Line2�� �����ϴ� 2���� ��ǥ��
	//   cross : �������� ��ǥ
	//   return: 0:����  1:����
	static int mathLineLineCross2D(double line1[2][2], double line2[2][2], double cross[2]);
	// ������ ������ ���� ���ϱ�
	//   line1 : Line�� �����ϴ� 2���� ��ǥ�� (line1[1][0] 2��° ��ǥ�� x��, line1[1][1] 2��° ��ǥ�� y��)
	//   sttP  : ������ �����ϴ� ��ǥ��
	//   endP  : ������ �����ϴ� ��ǥ��
	//   cross : �������� ��ǥ
	//   return: 0:�������� ���г��� ������  1:�������� ���г��� ������
	static int mathLineSegCross2D(double line[2][2], double sttP[2], double endP[2], double cross[2]);

	// 2���� polyline�� dOffset�Ÿ���ŭ ������ �ٱ���, �������� ũ�⸦ ������
	//   dOffset : Offset�Ÿ�(-:������ �������� Offset,  +:������ �ٱ������� Offset)
	//   nData   : polyLine�� ��ǥ ����
	// < polyLine: �Էµ� polyLine�� Offset���Ѽ� �Ѱ���
	//   nRotType: polyLine�� ȸ������(0:Auto,  1:�ݽð����,  2:�ð����)
	static bool mathOffsetOfPolyline(double dOffset, const int nData, double polyLine[][2], int nRotType = 0);

	// 3������ ���� ���� ������ ���
	//   P1, P2, P3 : ȣ�� �����ϴ� 3��(1,2,3����)
	// < CenterP    : ���� �߽���
	// < Radius     : ������
	//   RETURN     : -1:3���� ��� ���� ���, 0:3���� �����ΰ��, 1:��갡���� ���
	static int  mathCircleForm3Point(double p1[2], double p2[2], double p3[2], double CenterP[2], double& dRadius);
	// 3������ ���� ȣ�� ������ ���
	//   P1, P2, P3 : ȣ�� �����ϴ� 3��(1,2,3����)
	// < CenterP    : ���� �߽���
	// < Radius     : ������
	// < StartAngle : ���۰���
	// < SweepAngle : ���ΰ���(�ݽð���� +)
	//   RETURN     : -1:3���� ��� ���� ���, 0:3���� �����ΰ��, 1:��갡���� ���
	static int  mathArcForm3Point(double p1[2], double p2[2], double p3[2], double CenterP[2], double& dRadius, double& dStartAngle, double& dSweepAngle);
	static int  mathArcForm3Point(double p1[2], double p2[2], double p3[2], double CenterP[2], double& dRadius, double& dStartAngle, double& dSweepAngle, double& dBulge);

	static bool   GetCenter_Bulge(double line[2][2], double bulge, double centerP[2]);
	static double GetBulge(double line[2][2], double radius, bool inside);
	static int    GetRotationType(const int nData, double polyLine[][2]);//polyLine�� ȸ�������� �Ѱ���(1:�ݽð����,  2:�ð����)

	// �ش������ �������� ������� ���� ������ ���(����� ����̾����)�� ���� �� �߽����� ����Ͽ� �Ѱ���
	// < pointCenV, dVolm : ���� �߽����� ����
	//   nData, PontUnitList : ��� ����� 3D��ǥ �� ���� ����
	//   nAxis : ����� ����� ���� (1:xy���, 2:yz���, 3:xz���)
	static bool Get_VolmCenter(double pointCenV[3], double& dVolm, const int nData, double PontUnitList[][3], int nAxis = 1);

	// �� Polyline�� ���� ��ø���� üũ(���� �ǰų� �������� ��ĥ��� true, ��ġ�� ���� ���� ��� false) - Test�� ��뺸��(by SHIN 2010.05.31)
	static bool IsOverLap_Polyline(const int nData1, double polyLine1[][2], const int nData2, double polyLine2[][2]);

	// �������� �������� ����
	static bool IntersectTriangle(double RayOrigin[3], double RayDirection[3], double V0[3], double V1[3], double V2[3], double cross[3]);
	static bool IntersectTriangle(const dgn::math::CVector3& RayOrigin, const dgn::math::CVector3& RayDirection,
		const dgn::math::CVector3& V0, const dgn::math::CVector3& V1, const dgn::math::CVector3& V2, dgn::math::CVector3& cross);

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
	static int RealToInt(const double dVal, double dTol = 1e-6);

public:
	static double mathMomentOfInertiaRectangle(double dB, double dH, double de);
	static double mathMomentOfInertiaTriangle(double dB, double dH, double de);
	static double mathMomentOfInertiaCircle(double dr);

private:
	static bool IsNormalZero(const double& dVal);
};
#include "HeaderPost.h"
