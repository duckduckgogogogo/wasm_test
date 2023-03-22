#include "MathFunc.h"
#include <float.h>

#define __max(a,b) (((a) > (b)) ? (a) : (b))
#define __min(a,b) (((a) < (b)) ? (a) : (b))

using namespace dgn::math;
double CDgnMathFunc::m_NormalZero = 1.e-15;
double CDgnMathFunc::m_pi         = 4.0 * atan(1.0);
double CDgnMathFunc::m_trrad      = m_pi / 180.0;
double CDgnMathFunc::m_trang      = 180.0 / m_pi;
double CDgnMathFunc::m_VertTol    = 0.01;
double CDgnMathFunc::m_AngleZero  = 1.e-8; //degree

double CDgnMathFunc::GetAngle(double x1, double y1, double x2, double y2, bool bPosAngle)
{
	if ( CDgnMathFunc::mathLength(x1, y1, x2, y2) < 1e-3 )
	{
		return 0.0;
	}

	double dx = x1 - x2;
	double dy = y1 - y2;
	if ( fabs(dx) < CDgnMathFunc::Get_NormalZero() )
	{
		if ( dy >= 0.0 - CDgnMathFunc::Get_NormalZero() )
		{
			if ( bPosAngle )
			{
				return 270.0;
			}
			else
			{
				return -90.0;
			}
		}
		else
		{
			return 90.0;
		}
	}
	else if ( fabs(dy) < CDgnMathFunc::Get_NormalZero() )
	{
		if ( dx >= 0.0 - CDgnMathFunc::Get_NormalZero() )
		{
			if ( bPosAngle )
			{
				return 180.0;
			}
			else
			{
				return 0.0;
			}
		}
		else
		{
			return 0.0;
		}
	}
	else
	{
		double Angle = atan(dy / dx)*CDgnMathFunc::Get_trang();
		if ( !bPosAngle )
		{
			return Angle;
		}

		if ( dx < 0.0 && dy < 0.0 )
		{
			return fabs(Angle);
		}
		else if ( dx < 0.0 && dy > 0.0 )
		{
			return 360.0 - fabs(Angle);
		}
		else if ( dy > 0 && dx > 0 )
		{
			return 180.0 + fabs(Angle);
		}
		else
		{
			return 180.0 - fabs(Angle);
		}
	}
}

double CDgnMathFunc::CalcArea(const double& dDia)
{
	return Get_pi() * pow(dDia, 2.0) / 4.0;
}

double CDgnMathFunc::Get_pi()
{
	return m_pi;
}

double CDgnMathFunc::Get_trrad()
{
	return m_trrad;
}

double CDgnMathFunc::Get_trang()
{
	return m_trang;
}

double CDgnMathFunc::Get_NormalZero()
{
	return m_NormalZero;
}

bool CDgnMathFunc::IsNormalZero(const double& dVal)
{
	return ( fabs(dVal) < m_NormalZero ) ? true : false;
}

// 용도: counter clock wise인지 판단 1:CCW, -1:CW, 0:수평
// 담당: 선종복
int CDgnMathFunc::math_ccw(double ax, double ay, double bx, double by, double cx, double cy)
{
	double l =
		bx * cy - ay * bx - ax * cy -
		by * cx + ax * by + ay * cx;
	if ( fabs(l) < 1.0e-10 )
	{
		return 0;
	}
	else if ( l > 0.0 )
	{
		return 1;
	}
	else
	{
		return -1;
	}
}

void CDgnMathFunc::mathMxIdentity(mf_matrix m)
{
	static mf_matrix mtemp =
	{
		1.0, 0.0, 0.0, 0.0,
		0.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0, 0.0,
		0.0, 0.0, 0.0, 1.0
	};
	memcpy(m, mtemp, sizeof(mtemp));
}

void CDgnMathFunc::mathMxMult(mf_matrix m1, mf_matrix m2, mf_matrix mout)
{
	mf_matrix m;
	memset(m, 0, sizeof(m));

	for ( int i = 0; i < 4; i++ )
	{
		for ( int j = 0; j < 4; j++ )
		{
			for ( int k = 0; k < 4; k++ )
			{
				m[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}
	memcpy(mout, m, sizeof(m));
}

void CDgnMathFunc::mathMxCopy(mf_matrix min, mf_matrix mout)
{
	memcpy(mout, min, sizeof(mout));
}

void CDgnMathFunc::mathVecMult(mf_vector v, mf_matrix m, mf_vector vout)
{
	mf_vector vtemp;

	for ( int i = 0; i < 4; i++ )
	{
		vtemp[i] = 0;
		for ( int j = 0; j < 4; j++ )
		{
			vtemp[i] += v[j] * m[j][i];
		}
	}
	for ( int i = 0; i < 4; i++ )
	{
		vout[i] = vtemp[i];
	}
}

double CDgnMathFunc::mathVecDot(mf_vector v1, mf_vector v2)
{
	return ( v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2] );
}

void CDgnMathFunc::mathVecCross(mf_vector v1, mf_vector v2, mf_vector vout)
{
	mf_vector vtemp;
	vtemp[0] = v1[1] * v2[2] - v1[2] * v2[1];
	vtemp[1] = v1[2] * v2[0] - v1[0] * v2[2];
	vtemp[2] = v1[0] * v2[1] - v1[1] * v2[0];
	vtemp[3] = 1.0;
	memcpy(vout, vtemp, sizeof(vout));
}

void CDgnMathFunc::mathRotateX(double angle, double& rx, double& ry, double& rz)
{
	mf_matrix m;
	mf_vector v;
	double rad = angle * m_trrad;
	mathMxIdentity(m);
	double costh = mathCos(rad);
	double sinth = mathSin(rad);
	m[1][1] = costh;
	m[1][2] = sinth;
	m[2][1] = -sinth;
	m[2][2] = costh;
	v[0] = rx;
	v[1] = ry;
	v[2] = rz;
	v[3] = 1.0;
	mathVecMult(v, m, v);
	rx = v[0];
	ry = v[1];
	rz = v[2];
	if ( fabs(rx) < 1.e-10 ) rx = 0.0;
	if ( fabs(ry) < 1.e-10 ) ry = 0.0;
	if ( fabs(rz) < 1.e-10 ) rz = 0.0;
}

void CDgnMathFunc::mathRotateY(double angle, double& rx, double& ry, double& rz)
{
	mf_matrix m;
	mf_vector v;
	double rad = angle * m_trrad;
	mathMxIdentity(m);
	double costh = mathCos(rad);
	double sinth = mathSin(rad);
	m[0][0] = costh;
	m[0][2] = -sinth;
	m[2][0] = sinth;
	m[2][2] = costh;
	v[0] = rx;
	v[1] = ry;
	v[2] = rz;
	v[3] = 1.0;
	mathVecMult(v, m, v);
	rx = v[0];
	ry = v[1];
	rz = v[2];
	if ( fabs(rx) < 1.e-10 ) rx = 0.0;
	if ( fabs(ry) < 1.e-10 ) ry = 0.0;
	if ( fabs(rz) < 1.e-10 ) rz = 0.0;
}

void CDgnMathFunc::mathRotateZ(double angle, double& rx, double& ry, double& rz)
{
	mf_matrix m;
	mf_vector v;
	double rad = angle * m_trrad;
	mathMxIdentity(m);
	double costh = mathCos(rad);
	double sinth = mathSin(rad);
	m[0][0] = costh;
	m[0][1] = sinth;
	m[1][0] = -sinth;
	m[1][1] = costh;
	v[0] = rx;
	v[1] = ry;
	v[2] = rz;
	v[3] = 1.0;
	mathVecMult(v, m, v);
	rx = v[0];
	ry = v[1];
	rz = v[2];
	if ( fabs(rx) < 1.e-10 ) rx = 0.0;
	if ( fabs(ry) < 1.e-10 ) ry = 0.0;
	if ( fabs(rz) < 1.e-10 ) rz = 0.0;
}

void CDgnMathFunc::mathVec2Ang(double ux, double uy, double uz, double& rx, double& ry)
{
	rx = 0.0, ry = 0.0;
	double d1 = sqrt(uy*uy + uz*uz);
	if ( d1 >= 1.e-15 )
	{
		rx = asin(uy / d1)*m_trang;
		if ( uz < 0.0 )
		{
			rx = 180.0 - rx;
		}
	}
	double d2 = sqrt(ux*ux + uy*uy + uz*uz);
	if ( d2 >= 1.e-15 )
	{
		ry = asin(ux / d2)*m_trang;
	}
	rx = -rx;
}

void CDgnMathFunc::mathRotate(double angle, double ux, double uy, double uz, double& rx, double& ry, double& rz)
{
	if ( angle == 0.0 ) return;

	double d1 = sqrt(uy*uy + uz*uz);
	double phi = 0.0;
	if ( d1 >= 1.e-15 )
	{
		phi = asin(uy / d1)*m_trang;
		if ( uz < 0.0 )
		{
			phi = 180.0 - phi;
		}
	}

	double d2 = sqrt(ux*ux + uy * uy + uz * uz);
	double theta = 0.0;
	if ( d2 >= 1.e-15 )
	{
		theta = asin(ux / d2)*m_trang;
	}

	if ( phi   != 0.0 ) mathRotateX(phi, rx, ry, rz);
	if ( theta != 0.0 ) mathRotateY(-theta, rx, ry, rz);
	mathRotateZ(angle, rx, ry, rz);
	if ( theta != 0.0 ) mathRotateY(theta, rx, ry, rz);
	if ( phi   != 0.0 ) mathRotateX(-phi, rx, ry, rz);

	if ( fabs(rx) < 1.e-10 ) rx = 0.0;
	if ( fabs(ry) < 1.e-10 ) ry = 0.0;
	if ( fabs(rz) < 1.e-10 ) rz = 0.0;
}

void CDgnMathFunc::mathRotate(double angle, double px, double py, double pz, double ux, double uy, double uz, double& rx, double& ry, double& rz)
{
	if ( angle == 0.0 ) return;
	double rxx = rx - px;
	double ryy = ry - py;
	double rzz = rz - pz;
	mathRotate(angle, ux, uy, uz, rxx, ryy, rzz);
	rx = rxx + px;
	ry = ryy + py;
	rz = rzz + pz;
	if ( fabs(rx) < 1.e-10 ) rx = 0.0;
	if ( fabs(ry) < 1.e-10 ) ry = 0.0;
	if ( fabs(rz) < 1.e-10 ) rz = 0.0;
	return;
}

double CDgnMathFunc::mathDot(double xyz1[3], double xyz2[3], bool bNormal)
{
	double dblDot = xyz1[0] * xyz2[0] + xyz1[1] * xyz2[1] + xyz1[2] * xyz2[2];
	if ( bNormal )
	{
		if ( dblDot > 1.0 ) dblDot = 1.0;
		if ( dblDot < -1.0 ) dblDot = -1.0;
	}
	return dblDot;
}

double CDgnMathFunc::mathDot(double u1x, double u1y, double u2x, double u2y, bool bNormal)
{
	double dblDot = u1x * u2x + u1y * u2y;
	if ( bNormal )
	{
		if ( dblDot > 1.0 ) dblDot = 1.0;
		if ( dblDot < -1.0 ) dblDot = -1.0;
	}
	return dblDot;
}

void CDgnMathFunc::mathCross(double xyz1[3], double xyz2[3], double xyzout[3])
{
	double xyztemp[3];
	xyztemp[0] = xyz1[1] * xyz2[2] - xyz1[2] * xyz2[1];
	xyztemp[1] = xyz1[2] * xyz2[0] - xyz1[0] * xyz2[2];
	xyztemp[2] = xyz1[0] * xyz2[1] - xyz1[1] * xyz2[0];
	xyzout[0] = xyztemp[0];
	xyzout[1] = xyztemp[1];
	xyzout[2] = xyztemp[2];
}

// 내부에서 Normalize
double CDgnMathFunc::mathCrossAngleNormalize(double Vector1[3], double Vector2[3])
{
	double d1 = Vector1[0] * Vector1[0] + Vector1[1] * Vector1[1] + Vector1[2] * Vector1[2];
	double d2 = Vector2[0] * Vector2[0] + Vector2[1] * Vector2[1] + Vector2[2] * Vector2[2];
	d1 = sqrt(d1);
	d2 = sqrt(d2);
	double dTol = 1.0e-10;
	if ( fabs(d1) < dTol ) { /*ASSERT(0);*/ return 0; }
	if ( fabs(d2) < dTol ) { /*ASSERT(0);*/ return 0; }
	double NV1[3], NV2[3];
	for ( int i = 0; i < 3; i++ )
	{
		NV1[i] = Vector1[i] / d1;
		NV2[i] = Vector2[i] / d2;
	}
	return mathCrossAngle(NV1, NV2);
}

// Normalize 시키고 부호도 판별
// V1에서 V2로 반시계방향이면 +, 시계 방향이면 -
double CDgnMathFunc::mathCrossAngleNormalizeWithSign(double Vector1[3], double Vector2[3])
{
	double dAngle = mathCrossAngleNormalize(Vector1, Vector2);
	double VectorN[3];

	// 두선의 법선 구함
	// 법선의 z가 +이면 CCW, -이면 CW
	mathNormal(Vector1, Vector2, VectorN);

	double dTol = 1.0e-10;
	if ( fabs(VectorN[2]) <= dTol ) // Normal Vector가 XY 평면인 경우
	{
		if ( fabs(VectorN[0]) <= dTol )  // X 가 0이면 Y 좌표 본다.
		{
			if ( VectorN[1] >= 0.0 ) return dAngle;
			else return -dAngle;
		}
		else  // Y가 0이면 X좌표 본다.
		{
			if ( VectorN[0] >= 0.0 ) return dAngle;
			else return -dAngle;
		}
	}
	if ( VectorN[2] >= 0.0 )  // CCW
	{
		return dAngle;
	}
	return -dAngle;
}

double CDgnMathFunc::mathCrossAngle(double Vector1[3], double Vector2[3])
{// Normalized 된 벡터만 써야한다.
	double v = mathDot(Vector1, Vector2);
	double r = mathAcos(v);
	return r * m_trang;
}

double CDgnMathFunc::mathCrossAngle2D(double u1x, double u1y, double u2x, double u2y)
{
	double v = mathDot(u1x, u1y, u2x, u2y);
	double r = mathAcos(v);
	return r * m_trang;
}

double CDgnMathFunc::mathCrossAngle2D(double Vector1[3], double Vector2[3])
{// Normalized 된 2차원 벡터(z=0)만 써야한다.
	// ASSERT(Vector1[2] == 0 && Vector2[2] == 0);

	double v = mathDot(Vector1, Vector2);
	double r = mathAcos(v);
	double n[3];
	mathCross(Vector1, Vector2, n); // v1->v2의 normal Vector를 구한다.
	double np[3] = { 0,0,1 };
	if ( mathDot(n, np) > 0 )
	{
		return r * m_trang;
	}
	else
	{
		return -r * m_trang;
	}
}

double CDgnMathFunc::mathLength(double dx, double dy, double dz)
{
	return sqrt(dx*dx + dy * dy + dz * dz);
}

double CDgnMathFunc::mathLength(double dxi, double dyi, double dxj, double dyj)
{
	double dblLength = ( dxi - dxj )*( dxi - dxj ) + ( dyi - dyj )*( dyi - dyj );
	if ( dblLength < 0.0 ) dblLength = 0.0;
	return sqrt(dblLength);
}

double CDgnMathFunc::mathLength(double dxi, double dyi, double dzi, double dxj, double dyj, double dzj)
{
	double dblLength = ( dxi - dxj )*( dxi - dxj ) + ( dyi - dyj )*( dyi - dyj ) + ( dzi - dzj )*( dzi - dzj );
	if ( dblLength < 0.0 ) dblLength = 0.0;
	return sqrt(dblLength);
}

double CDgnMathFunc::mathLength(double xyz1[3], double xyz2[3])
{
	return CDgnMathFunc::mathLength(xyz1[0], xyz1[1], xyz1[2], xyz2[0], xyz2[1], xyz2[2]);
}

double CDgnMathFunc::mathArea(double xyz1[3], double xyz2[3], double xyz3[3])
{
	double Length1 = mathLength(xyz1[0] - xyz2[0], xyz1[1] - xyz2[1], xyz1[2] - xyz2[2]);
	double Length2 = mathLength(xyz2[0] - xyz3[0], xyz2[1] - xyz3[1], xyz2[2] - xyz3[2]);
	double Length3 = mathLength(xyz3[0] - xyz1[0], xyz3[1] - xyz1[1], xyz3[2] - xyz1[2]);
	double s  = ( Length1 + Length2 + Length3 ) / 2.0;
	double ss = s * ( s - Length1 )*( s - Length2 )*( s - Length3 );
	if ( ss < 0.0 ) ss = 0.0; // JBSEON-20060313 2006-03-13 메일 (from 이대근 : kabuto-rev-zero.mcb)
	return sqrt(ss);
}

double CDgnMathFunc::mathArea2d(double x1, double y1, double x2, double y2, double x3, double y3)
{
	double Length1 = mathLength(x1, y1, x2, y2);
	double Length2 = mathLength(x2, y2, x3, y3);
	double Length3 = mathLength(x3, y3, x1, y1);
	double s  = ( Length1 + Length2 + Length3 ) / 2.0;
	double ss = s * ( s - Length1 )*( s - Length2 )*( s - Length3 );
	if ( ss < 0.0 ) ss = 0.0; // JBSEON-20060313 2006-03-13 메일 (from 이대근 : kabuto-rev-zero.mcb)
	return sqrt(ss);
}

double CDgnMathFunc::mathArea2d(int nPoint, double x[], double y[])
{
	double dArea = 0.0;
	double dx, dy, ymin;
	for ( int i = 0; i < nPoint; i++ )
	{
		if ( i == 0 ) ymin = y[i];
		if ( y[i] < ymin ) ymin = y[i];
	}
	for ( int i = 0; i < nPoint; i++ )
	{
		int j = ( i + 1 ) % nPoint;
		dx = x[i] - x[j];
		dy = ( y[i] + y[j] - 2 * ymin )*0.5;
		dArea = dArea + dx * dy;
	}
	return fabs(dArea);
}

double CDgnMathFunc::mathVolume(int nVertex, double xyz[][3])
{
	switch ( nVertex )
	{
	case 4: return Solid_4(xyz);
	case 6: return Solid_6(xyz);
	case 8: return Solid_8(xyz);
	}
	return 0.0;
}

double CDgnMathFunc::Solid_4(double xyz[][3])
{
	const int NumNode = 4;
	double T_DET_J = 0.0, DET_J;
	double DX_S = 0.0, DY_S = 0.0, DZ_S = 0.0;
	double DX_T = 0.0, DY_T = 0.0, DZ_T = 0.0;
	double DX_U = 0.0, DY_U = 0.0, DZ_U = 0.0;
	double DN_S[NumNode] = { -1.0, 1.0, 0.0, 0.0 };
	double DN_T[NumNode] = { -1.0, 0.0, 1.0, 0.0 };
	double DN_U[NumNode] = { -1.0, 0.0, 0.0, 1.0 };

	for ( int i = 0; i < NumNode; i++ )
	{
		DX_S += DN_S[i] * xyz[i][0];
		DY_S += DN_S[i] * xyz[i][1];
		DZ_S += DN_S[i] * xyz[i][2];
		DX_T += DN_T[i] * xyz[i][0];
		DY_T += DN_T[i] * xyz[i][1];
		DZ_T += DN_T[i] * xyz[i][2];
		DX_U += DN_U[i] * xyz[i][0];
		DY_U += DN_U[i] * xyz[i][1];
		DZ_U += DN_U[i] * xyz[i][2];
	}
	DET_J = DX_S * ( DY_T * DZ_U - DZ_T * DY_U )
		- DY_S * ( DX_T * DZ_U - DZ_T * DX_U )
		+ DZ_S * ( DX_T * DY_U - DY_T * DX_U );
	T_DET_J = T_DET_J + DET_J / 6.0;

	return T_DET_J;
}

double CDgnMathFunc::Solid_6(double xyz[][3])
{
	const int NumNode = 6;
	double T_DET_J = 0.0, DET_J;
	double  SS[8] = { +0.33333333, +0.20000000, +0.60000000, +0.20000000, +0.33333333, +0.20000000, +0.60000000, +0.20000000 };
	double  TT[8] = { +0.33333333, +0.20000000, +0.20000000, +0.60000000, +0.33333333, +0.20000000, +0.20000000, +0.60000000 };
	double  UU[8] = { -0.57735027, -0.57735027, -0.57735027, -0.57735027, +0.57735027, +0.57735027, +0.57735027, +0.57735027 };
	double WHT[8] = { -0.28125000, +0.26041667, +0.26041667, +0.26041667, -0.28125000, +0.26041667, +0.26041667, +0.26041667 };
	double DX_S, DY_S, DZ_S;
	double DX_T, DY_T, DZ_T;
	double DX_U, DY_U, DZ_U;
	double DN_S[NumNode], DN_T[NumNode], DN_U[NumNode];
	double SV, TV, UV;

	for ( int i = 0; i < 8; i++ )
	{
		SV = SS[i];
		TV = TT[i];
		UV = UU[i];
		DX_S = 0.0;
		DY_S = 0.0;
		DZ_S = 0.0;
		DX_T = 0.0;
		DY_T = 0.0;
		DZ_T = 0.0;
		DX_U = 0.0;
		DY_U = 0.0;
		DZ_U = 0.0;
		DN_S[0] = -( 1.0 / 2.0 )*( 1.0 - UV );
		DN_S[1] = +( 1.0 / 2.0 )*( 1.0 - UV );
		DN_S[2] = 0.0;
		DN_S[3] = -( 1.0 / 2.0 )*( 1.0 + UV );
		DN_S[4] = +( 1.0 / 2.0 )*( 1.0 + UV );
		DN_S[5] = 0.0;
		DN_T[0] = -( 1.0 / 2.0 )*( 1.0 - UV );
		DN_T[1] = 0.0;
		DN_T[2] = +( 1.0 / 2.0 )*( 1.0 - UV );
		DN_T[3] = -( 1.0 / 2.0 )*( 1.0 + UV );
		DN_T[4] = 0.0;
		DN_T[5] = +( 1.0 / 2.0 )*( 1.0 + UV );
		DN_U[0] = -( 1.0 / 2.0 )*( 1.0 - SV - TV );
		DN_U[1] = -( 1.0 / 2.0 )*( SV );
		DN_U[2] = -( 1.0 / 2.0 )*( TV );
		DN_U[3] = +( 1.0 / 2.0 )*( 1.0 - SV - TV );
		DN_U[4] = +( 1.0 / 2.0 )*( SV );
		DN_U[5] = +( 1.0 / 2.0 )*( TV );
		for ( int j = 0; j < NumNode; j++ )
		{
			DX_S += DN_S[j] * xyz[j][0];
			DY_S += DN_S[j] * xyz[j][1];
			DZ_S += DN_S[j] * xyz[j][2];
			DX_T += DN_T[j] * xyz[j][0];
			DY_T += DN_T[j] * xyz[j][1];
			DZ_T += DN_T[j] * xyz[j][2];
			DX_U += DN_U[j] * xyz[j][0];
			DY_U += DN_U[j] * xyz[j][1];
			DZ_U += DN_U[j] * xyz[j][2];
		}
		DET_J = DX_S * ( DY_T * DZ_U - DZ_T * DY_U )
			- DY_S * ( DX_T * DZ_U - DZ_T * DX_U )
			+ DZ_S * ( DX_T * DY_U - DY_T * DX_U );
		T_DET_J = T_DET_J + DET_J * WHT[i];
	}
	return T_DET_J;
}

double CDgnMathFunc::Solid_8(double xyz[][3])
{
	const int NumNode = 8;
	double T_DET_J = 0.0, DET_J;
	double  SS[8] = { -0.5773503, +0.5773503, +0.5773503, -0.5773503, -0.5773503, +0.5773503, +0.5773503, -0.5773503 };
	double  TT[8] = { -0.5773503, -0.5773503, +0.5773503, +0.5773503, -0.5773503, -0.5773503, +0.5773503, +0.5773503 };
	double  UU[8] = { -0.5773503, -0.5773503, -0.5773503, -0.5773503, +0.5773503, +0.5773503, +0.5773503, +0.5773503 };
	double DX_S, DY_S, DZ_S;
	double DX_T, DY_T, DZ_T;
	double DX_U, DY_U, DZ_U;
	double DN_S[NumNode], DN_T[NumNode], DN_U[NumNode];
	double SV, TV, UV, SM, SP, TM, TP, UM, UP;

	for ( int i = 0; i < 8; i++ )
	{
		SV = SS[i];
		TV = TT[i];
		UV = UU[i];
		DX_S = 0.0;
		DY_S = 0.0;
		DZ_S = 0.0;
		DX_T = 0.0;
		DY_T = 0.0;
		DZ_T = 0.0;
		DX_U = 0.0;
		DY_U = 0.0;
		DZ_U = 0.0;
		SM = ( 1.0 - SV );
		SP = ( 1.0 + SV );
		TM = ( 1.0 - TV );
		TP = ( 1.0 + TV );
		UM = ( 1.0 - UV );
		UP = ( 1.0 + UV );
		DN_S[0] = -( TM*UM ) / 8.0;
		DN_S[1] = +( TM*UM ) / 8.0;
		DN_S[2] = +( TP*UM ) / 8.0;
		DN_S[3] = -( TP*UM ) / 8.0;
		DN_S[4] = -( TM*UP ) / 8.0;
		DN_S[5] = +( TM*UP ) / 8.0;
		DN_S[6] = +( TP*UP ) / 8.0;
		DN_S[7] = -( TP*UP ) / 8.0;
		DN_T[0] = -( SM*UM ) / 8.0;
		DN_T[1] = -( SP*UM ) / 8.0;
		DN_T[2] = +( SP*UM ) / 8.0;
		DN_T[3] = +( SM*UM ) / 8.0;
		DN_T[4] = -( SM*UP ) / 8.0;
		DN_T[5] = -( SP*UP ) / 8.0;
		DN_T[6] = +( SP*UP ) / 8.0;
		DN_T[7] = +( SM*UP ) / 8.0;
		DN_U[0] = -( SM*TM ) / 8.0;
		DN_U[1] = -( SP*TM ) / 8.0;
		DN_U[2] = -( SP*TP ) / 8.0;
		DN_U[3] = -( SM*TP ) / 8.0;
		DN_U[4] = +( SM*TM ) / 8.0;
		DN_U[5] = +( SP*TM ) / 8.0;
		DN_U[6] = +( SP*TP ) / 8.0;
		DN_U[7] = +( SM*TP ) / 8.0;
		for ( int j = 0; j < NumNode; j++ )
		{
			DX_S += DN_S[j] * xyz[j][0];
			DY_S += DN_S[j] * xyz[j][1];
			DZ_S += DN_S[j] * xyz[j][2];
			DX_T += DN_T[j] * xyz[j][0];
			DY_T += DN_T[j] * xyz[j][1];
			DZ_T += DN_T[j] * xyz[j][2];
			DX_U += DN_U[j] * xyz[j][0];
			DY_U += DN_U[j] * xyz[j][1];
			DZ_U += DN_U[j] * xyz[j][2];
		}
		DET_J = DX_S * ( DY_T * DZ_U - DZ_T * DY_U )
			- DY_S * ( DX_T * DZ_U - DZ_T * DX_U )
			+ DZ_S * ( DX_T * DY_U - DY_T * DX_U );
		T_DET_J += DET_J;
	}
	return T_DET_J;
}

void CDgnMathFunc::mathNormal(double Vector1[3], double Vector2[3], double VectorN[3])
{
	double Vector[3];
	Vector[0] = Vector1[1] * Vector2[2] - Vector1[2] * Vector2[1];
	Vector[1] = Vector1[2] * Vector2[0] - Vector1[0] * Vector2[2];
	Vector[2] = Vector1[0] * Vector2[1] - Vector1[1] * Vector2[0];
	VectorN[0] = Vector[0];
	VectorN[1] = Vector[1];
	VectorN[2] = Vector[2];
}

void CDgnMathFunc::mathNormal(double Vector1[3], double Vector2[3], double Vector3[3], double VectorN[3])
{
	double Vector[3];
	Vector[0] = ( Vector2[1] - Vector1[1] )*( Vector3[2] - Vector1[2] ) - ( Vector2[2] - Vector1[2] )*( Vector3[1] - Vector1[1] );
	Vector[1] = ( Vector2[2] - Vector1[2] )*( Vector3[0] - Vector1[0] ) - ( Vector2[0] - Vector1[0] )*( Vector3[2] - Vector1[2] );
	Vector[2] = ( Vector2[0] - Vector1[0] )*( Vector3[1] - Vector1[1] ) - ( Vector2[1] - Vector1[1] )*( Vector3[0] - Vector1[0] );
	VectorN[0] = Vector[0];
	VectorN[1] = Vector[1];
	VectorN[2] = Vector[2];
}

bool CDgnMathFunc::mathNormalize(double xyz[3], double xyzn[3])
{
	double Length = mathLength(xyz[0], xyz[1], xyz[2]);
	if ( IsNormalZero(Length) ) return false;
	xyzn[0] = xyz[0] / Length;
	xyzn[1] = xyz[1] / Length;
	xyzn[2] = xyz[2] / Length;
	return true;
}

bool CDgnMathFunc::mathNormalize(double dx, double dy, double dz, double& dxn, double& dyn, double& dzn)
{
	double Length = mathLength(dx, dy, dz);
	dxn = 0.0;
	dyn = 0.0;
	dzn = 0.0;
	if ( IsNormalZero(Length) ) return false;
	dxn = dx / Length;
	dyn = dy / Length;
	dzn = dz / Length;
	return true;
}

bool CDgnMathFunc::mathNormalize(double dx, double dy, double& dxn, double& dyn)
{
	double Length = mathLength(dx, dy);
	dxn = 0.0;
	dyn = 0.0;
	if ( IsNormalZero(Length) ) return false;
	dxn = dx / Length;
	dyn = dy / Length;
	return true;
}

void CDgnMathFunc::mathPLCrossPoint(double a, double b, double c, double l, double m, double n,
	double x, double y, double z, double& rpx, double& rpy, double& rpz)
{
	double t = ( ( x - a )*l + ( y - b )*m + ( z - c )*n ) / ( l*l + m * m + n * n );
	rpx = a + t * l;
	rpy = b + t * m;
	rpz = c + t * n;
}

void CDgnMathFunc::mathPLCrossPoint(double a, double b, double c, double l, double m, double n,
	double x, double y, double z, double& rpx, double& rpy, double& rpz, double& t)
{
	t = ( ( x - a )*l + ( y - b )*m + ( z - c )*n ) / ( l*l + m * m + n * n );
	rpx = a + t * l;
	rpy = b + t * m;
	rpz = c + t * n;
}

bool CDgnMathFunc::mathPlaneEquation(double p1[3], double p2[3], double p3[3], double& a, double& b, double& c, double & h)
{
	a = 0.0, b = 0.0, c = 0.0, h = 0.0;

	double aa, bb, cc, dd;
	aa = p1[1] * ( p2[2] - p3[2] ) - p2[1] * ( p1[2] - p3[2] ) + p3[1] * ( p1[2] - p2[2] );
	bb = p1[0] * ( p3[2] - p2[2] ) + p2[0] * ( p1[2] - p3[2] ) + p3[0] * ( p2[2] - p1[2] );
	cc = p1[0] * ( p2[1] - p3[1] ) - p2[0] * ( p1[1] - p3[1] ) + p3[0] * ( p1[1] - p2[1] );
	dd = CDgnMathFunc::mathSqrt(aa*aa + bb * bb + cc * cc);
	if ( IsNormalZero(dd) ) return false;
	a = aa / dd;
	b = bb / dd;
	c = cc / dd;
	h = -( a*p1[0] + b * p1[1] + c * p1[2] );
	return true;
}

double CDgnMathFunc::mathDistanceToPoint(double p1[3], double p2[3])
{
	double dx = p1[0] - p2[0];
	double dy = p1[1] - p2[1];
	double dz = p1[2] - p2[2];
	return sqrt(dx*dx + dy * dy + dz * dz);
}

double CDgnMathFunc::mathDistanceToLine(double line_i[3], double line_j[3], double point[3])
{
	double a, b, c, l, m, n, x, y, z;
	double dx, dy, dz, px, py, pz;

	x = point[0], y = point[1], z = point[2];
	a = line_i[0], b = line_i[1], c = line_i[2];
	dx = line_j[0] - a, dy = line_j[1] - b, dz = line_j[2] - c;

	double dLength = mathLength(dx, dy, dz);
	mathNormalize(dx, dy, dz, l, m, n);
	double t;
	if ( IsNormalZero(l) && IsNormalZero(m) && IsNormalZero(n) )
	{
		t = 0.0;
	}
	else
	{
		t = ( ( x - a )*l + ( y - b )*m + ( z - c )*n ) / ( l*l + m * m + n * n );
	}

	if ( t <= 0.0 )
	{
		px = line_i[0], py = line_i[1], pz = line_i[2];
	}
	else if ( t >= dLength )
	{
		px = line_j[0], py = line_j[1], pz = line_j[2];
	}
	else
	{
		px = a + t * l;
		py = b + t * m;
		pz = c + t * n;
	}
	return mathLength(point[0] - px, point[1] - py, point[2] - pz);
}

double CDgnMathFunc::mathDistanceToLine2D(double line_i_x, double line_i_y, double line_j_x, double line_j_y, double point_x, double point_y)
{
	double a, b, l, m, x, y;
	double dx, dy, px, py;

	x = point_x, y = point_y;
	a = line_i_x, b = line_i_y;
	dx = line_j_x - a, dy = line_j_y - b;

	double dLength = mathLength(dx, dy);
	mathNormalize(dx, dy, l, m);
	double t;
	if ( IsNormalZero(l) && IsNormalZero(m) )
	{
		t = 0.0;
	}
	else
	{
		t = ( ( x - a )*l + ( y - b )*m ) / ( l*l + m * m );
	}
	if ( t <= 0.0 )
	{
		px = line_i_x, py = line_i_y;
	}
	else if ( t >= dLength )
	{
		px = line_j_x, py = line_j_y;
	}
	else
	{
		px = a + t * l;
		py = b + t * m;
	}
	return mathLength(point_x - px, point_y - py);
}

double CDgnMathFunc::mathDistanceToPlane(double a, double b, double c, double h, double x, double y, double z)
{
	double aa = fabs(a*x + b * y + c * z + h);
	double dd = CDgnMathFunc::mathSqrt(a*a + b * b + c * c);
	if ( IsNormalZero(dd) ) return -1.0;
	return aa / dd;
}

double CDgnMathFunc::mathDistanceToPlaneWithSign(double a, double b, double c, double h, double x, double y, double z)
{
	double aa = a * x + b * y + c * z + h;
	double dd = CDgnMathFunc::mathSqrt(a*a + b * b + c * c);
	if ( IsNormalZero(dd) ) return -1.0;
	return aa / dd;
}

bool CDgnMathFunc::mathIncludePointInPlane(double p1[3], double p2[3], double p3[3], double point[3], double Tol)
{
	double a, b, c, h;
	double dDistance = 0.0;
	if ( !mathPlaneEquation(p1, p2, p3, a, b, c, h) ) return false;

	dDistance = mathDistanceToPlane(a, b, c, h, point[0], point[1], point[2]);
	if ( fabs(dDistance) < Tol )
	{
		return true;
	}
	else
	{
		return false;
	}
}

double CDgnMathFunc::mathDistanceFromIntersectPointToLine(double line_i[3], double line_j[3], double point[3])
{
	double pInt[3] = { 0.0, };
	mathIntersectPointToLine(line_i, line_j, point, pInt);
	return mathLength(point[0] - pInt[0], point[1] - pInt[1], point[2] - pInt[2]);
}

bool CDgnMathFunc::mathIntersectLine2D(double p1i_x, double p1i_y, double p1j_x, double p1j_y, double p2i_x, double p2i_y, double p2j_x, double p2j_y, double& pInts_x, double& pInts_y, bool bCheckOnLine)
{
	pInts_x = 0.0;
	pInts_y = 0.0;
	double u1, v1, u2, v2;
	if ( !mathNormalize(p1j_x - p1i_x, p1j_y - p1i_y, u1, v1) ) return false;
	if ( !mathNormalize(p2j_x - p2i_x, p2j_y - p2i_y, u2, v2) ) return false;
	double d = sqrt(( u1*v2 - u2 * v1 )*( u1*v2 - u2 * v1 ));

	if ( fabs(d) < 1.e-4 ) return false;

	double q = u2 * v1 - u1 * v2;
	double t;
	double eps = 1.e-7;
	if ( fabs(q) < eps ) return false;
	t = ( u1*( p2i_y - p1i_y ) - v1 * ( p2i_x - p1i_x ) ) / q;

	if ( bCheckOnLine )
	{
		double dLength = mathLength(p2i_x, p2i_y, p2j_x, p2j_y);
		if ( t<-dLength * eps || t> dLength + dLength * eps ) return false;
		dLength = mathLength(p1i_x, p1i_y, p1j_x, p1j_y);
		q = u1 * v2 - u2 * v1;
		double tt = ( u2*( p1i_y - p2i_y ) - v2 * ( p1i_x - p2i_x ) ) / q;
		if ( tt<-dLength * eps || tt> dLength + dLength * eps ) return false;
	}

	// pl2선상의 교점을 구한다.
	pInts_x = p2i_x + u2 * t;
	pInts_y = p2i_y + v2 * t;
	return true;
}

bool CDgnMathFunc::mathIntersectLinePMCv(double p1i_x, double p1i_y, double p1j_x, double p1j_y, double p2i_x, double p2i_y, double p2j_x, double p2j_y, double& pInts_x, double& pInts_y)
{
	pInts_x = 0.0;
	pInts_y = 0.0;
	double u1, v1, u2, v2;
	if ( !mathNormalize(p1j_x - p1i_x, p1j_y - p1i_y, u1, v1) ) return false;
	if ( !mathNormalize(p2j_x - p2i_x, p2j_y - p2i_y, u2, v2) ) return false;
	double d = sqrt(( u1*v2 - u2 * v1 )*( u1*v2 - u2 * v1 ));
	if ( fabs(d) < 1.e-7 ) return false;

	double q = u2 * v1 - u1 * v2;
	double eps = 1.e-7;
	if ( fabs(q) < eps ) return false;

	double t = ( u1*( p2i_y - p1i_y ) - v1 * ( p2i_x - p1i_x ) ) / q;
	pInts_x = p2i_x + u2 * t;
	pInts_y = p2i_y + v2 * t;
	return true;
}

bool CDgnMathFunc::mathIntersectPointToLine(double line_i[3], double line_j[3], double point[3], double rpInt[3])
{
	double x = point[0];
	double y = point[1];
	double z = point[2];

	double a = line_i[0];
	double b = line_i[1];
	double c = line_i[2];

	double dx = line_j[0] - a;
	double dy = line_j[1] - b;
	double dz = line_j[2] - c;

	double l = 0.0, m = 0.0, n = 0.0;
	mathNormalize(dx, dy, dz, l, m, n);

	double t = 0.0;

	if ( !(l*l + m * m + n * n == 0.0) )
	{
		t = ( ( x - a )*l + ( y - b )*m + ( z - c )*n ) / ( l*l + m * m + n * n );
	}

	if ( IsNormalZero(l) && IsNormalZero(m) && IsNormalZero(n) )
	{
		t = 0.0;
	}

	rpInt[0] = a + t * l;
	rpInt[1] = b + t * m;
	rpInt[2] = c + t * n;
	return true;
}

bool CDgnMathFunc::mathIntersectPointToPlane(double LineVector[3], double LinePoint[3], double PlaneNormal[3], double PlanePoint[3], double pInts[3])
{
	pInts[0] = 0.0, pInts[1] = 0.0, pInts[2] = 0.0;
	if ( !mathNormalize(LineVector, LineVector) ) return false;
	if ( !mathNormalize(PlaneNormal, PlaneNormal) ) return false;

	double v = LineVector[0] * PlaneNormal[0] + LineVector[1] * PlaneNormal[1] + LineVector[2] * PlaneNormal[2];
	double a =
		PlaneNormal[0] * PlanePoint[0] + PlaneNormal[1] * PlanePoint[1] + PlaneNormal[2] * PlanePoint[2] -
		PlaneNormal[0] * LinePoint[0]  - PlaneNormal[1] * LinePoint[1]  - PlaneNormal[2] * LinePoint[2];

	/****************************************************
	@@ 주어진 직선과 평면이 평행하다
	**/
	if ( fabs(v) <= 1.e-7 ) return false;
	double t = a / v;

	pInts[0] = LinePoint[0] + LineVector[0] * t;
	pInts[1] = LinePoint[1] + LineVector[1] * t;
	pInts[2] = LinePoint[2] + LineVector[2] * t;
	return true;
}

bool CDgnMathFunc::mathIntersectPointToPlane(double LinePoint1[3], double LinePoint2[3], double PlaneNormal[3], double PlanePoint[3], double pInts[3], double& dDistance)
{
	pInts[0] = 0.0, pInts[1] = 0.0, pInts[2] = 0.0;
	dDistance = 0.0;

	double LineVector[3];
	double dx = LinePoint2[0] - LinePoint1[0];
	double dy = LinePoint2[1] - LinePoint1[1];
	double dz = LinePoint2[2] - LinePoint1[2];
	double dLength = mathLength(dx, dy, dz);
	if ( dLength < 1.e-10 ) return false;

	LineVector[0] = dx / dLength, LineVector[1] = dy / dLength, LineVector[2] = dz / dLength;
	if ( !mathNormalize(PlaneNormal, PlaneNormal) ) return false;

	double v = LineVector[0] * PlaneNormal[0] + LineVector[1] * PlaneNormal[1] + LineVector[2] * PlaneNormal[2];
	double a =
		PlaneNormal[0] * PlanePoint[0] + PlaneNormal[1] * PlanePoint[1] + PlaneNormal[2] * PlanePoint[2]-
		PlaneNormal[0] * LinePoint1[0] - PlaneNormal[1] * LinePoint1[1] - PlaneNormal[2] * LinePoint1[2];

	/****************************************************
	@@ 주어진 직선과 평면이 평행하다
	**/
	if ( fabs(v) <= 1.e-7 ) return false;

	double t = a / v;
	pInts[0] = LinePoint1[0] + LineVector[0] * t;
	pInts[1] = LinePoint1[1] + LineVector[1] * t;
	pInts[2] = LinePoint1[2] + LineVector[2] * t;

	if ( t < 0.0 ) dDistance = fabs(t);
	if ( t > dLength ) dDistance = fabs(t) - dLength;
	return true;
}

double CDgnMathFunc::mathDistanceLine(double pl1_i[3], double pl1_j[3], double pl2_i[3], double pl2_j[3])
{
	double dblDistance = -1.0;
	double u1x = pl1_j[0] - pl1_i[0];
	double u1y = pl1_j[1] - pl1_i[1];
	double u1z = pl1_j[2] - pl1_i[2];
	double u2x = pl2_j[0] - pl2_i[0];
	double u2y = pl2_j[1] - pl2_i[1];
	double u2z = pl2_j[2] - pl2_i[2];
	double u1, v1, w1, u2, v2, w2;
	if ( !mathNormalize(u1x, u1y, u1z, u1, v1, w1) ) return dblDistance;
	if ( !mathNormalize(u2x, u2y, u2z, u2, v2, w2) ) return dblDistance;

	double d = sqrt(( v1*w2 - v2 * w1 )*( v1*w2 - v2 * w1 ) + ( w1*u2 - w2 * u1 )*( w1*u2 - w2 * u1 ) + ( u1*v2 - u2 * v1 )*( u1*v2 - u2 * v1 ));
	if ( fabs(d) < 1.e-7 ) return dblDistance;

	dblDistance = fabs(( ( pl2_i[0] - pl1_i[0] )*( v1*w2 - v2 * w1 ) - ( pl2_i[1] - pl1_i[1] )*( u1*w2 - u2 * w1 ) + ( pl2_i[2] - pl1_i[2] )*( u1*v2 - u2 * v1 ) ) / d);
	if ( dblDistance < 0.0 ) dblDistance = 0.0;
	return dblDistance;
}

bool CDgnMathFunc::mathIntersectLine(double pl1_i[3], double pl1_j[3], double pl2_i[3], double pl2_j[3], double Tolerance, double& dblDistance, double pInts[3])
{
	pInts[0] = 0.0;
	pInts[1] = 0.0;
	pInts[2] = 0.0;
	double u1x = pl1_j[0] - pl1_i[0];
	double u1y = pl1_j[1] - pl1_i[1];
	double u1z = pl1_j[2] - pl1_i[2];
	double u2x = pl2_j[0] - pl2_i[0];
	double u2y = pl2_j[1] - pl2_i[1];
	double u2z = pl2_j[2] - pl2_i[2];
	double u1, v1, w1, u2, v2, w2;
	if ( !mathNormalize(u1x, u1y, u1z, u1, v1, w1) ) return false;
	if ( !mathNormalize(u2x, u2y, u2z, u2, v2, w2) ) return false;
	double d = sqrt(( v1*w2 - v2 * w1 )*( v1*w2 - v2 * w1 ) + ( w1*u2 - w2 * u1 )*( w1*u2 - w2 * u1 ) + ( u1*v2 - u2 * v1 )*( u1*v2 - u2 * v1 ));

	if ( fabs(d) < 1.e-4 ) return false;
	dblDistance = fabs(( ( pl2_i[0] - pl1_i[0] )*( v1*w2 - v2 * w1 ) - ( pl2_i[1] - pl1_i[1] )*( u1*w2 - u2 * w1 ) + ( pl2_i[2] - pl1_i[2] )*( u1*v2 - u2 * v1 ) ) / d);
	if ( dblDistance > Tolerance ) return false;
	if ( dblDistance < 0.0 ) dblDistance = 0.0;

	double q1 = u2 * v1 - u1 * v2;
	double q2 = u2 * w1 - u1 * w2;
	double q3 = v2 * w1 - v1 * w2;
	double t;

	double eps = 1.e-7;
	double maxq = __max(__max(fabs(q1), fabs(q2)), fabs(q3));
	if ( maxq < eps ) return false;
	if ( maxq == fabs(q1) ) t = ( u1*( pl2_i[1] - pl1_i[1] ) - v1 * ( pl2_i[0] - pl1_i[0] ) ) / q1;
	else if ( maxq == fabs(q2) ) t = ( u1*( pl2_i[2] - pl1_i[2] ) - w1 * ( pl2_i[0] - pl1_i[0] ) ) / q2;
	else if ( maxq == fabs(q3) ) t = ( v1*( pl2_i[2] - pl1_i[2] ) - w1 * ( pl2_i[1] - pl1_i[1] ) ) / q3;
	else
	{
		return false;
	}

	// pl2선상의 교점을 구한다.
	double p2x = pl2_i[0] + u2 * t;
	double p2y = pl2_i[1] + v2 * t;
	double p2z = pl2_i[2] + w2 * t;

	pInts[0] = p2x;
	pInts[1] = p2y;
	pInts[2] = p2z;
	return true;
}

bool CDgnMathFunc::mathIntersectLine2(double pl1_i[3], double pl1_j[3], double pl2_i[3], double pl2_j[3], double Tolerance, double& dblDistance, double pInts[3])
{
	pInts[0] = 0.0;
	pInts[1] = 0.0;
	pInts[2] = 0.0;
	double u1x = pl1_j[0] - pl1_i[0];
	double u1y = pl1_j[1] - pl1_i[1];
	double u1z = pl1_j[2] - pl1_i[2];
	double u2x = pl2_j[0] - pl2_i[0];
	double u2y = pl2_j[1] - pl2_i[1];
	double u2z = pl2_j[2] - pl2_i[2];
	double u1, v1, w1, u2, v2, w2;
	double dL1Length = mathLength(u1x, u1y, u1z);
	double dL2Length = mathLength(u2x, u2y, u2z);

	if ( !mathNormalize(u1x, u1y, u1z, u1, v1, w1) ) return false;
	if ( !mathNormalize(u2x, u2y, u2z, u2, v2, w2) ) return false;
	double d = sqrt(( v1*w2 - v2 * w1 )*( v1*w2 - v2 * w1 ) + ( w1*u2 - w2 * u1 )*( w1*u2 - w2 * u1 ) + ( u1*v2 - u2 * v1 )*( u1*v2 - u2 * v1 ));

	if ( fabs(d) < 1.e-6 ) return false;
	dblDistance = fabs(( ( pl2_i[0] - pl1_i[0] )*( v1*w2 - v2 * w1 ) - ( pl2_i[1] - pl1_i[1] )*( u1*w2 - u2 * w1 ) + ( pl2_i[2] - pl1_i[2] )*( u1*v2 - u2 * v1 ) ) / d);
	if ( dblDistance > Tolerance ) return false;
	if ( dblDistance < 0.0 )dblDistance = 0.0;

	double q1 = u1 * v2 - u2 * v1;
	double q2 = u1 * w2 - u2 * w1;
	double q3 = v1 * w2 - v2 * w1;
	double t;
	double eps = 1.e-7;
	double maxq = __max(__max(fabs(q1), fabs(q2)), fabs(q3));
	if ( maxq < eps ) return false;
	if ( maxq == fabs(q1) ) t = ( u2*( pl1_i[1] - pl2_i[1] ) - v2 * ( pl1_i[0] - pl2_i[0] ) ) / q1;
	else if ( maxq == fabs(q2) ) t = ( u2*( pl1_i[2] - pl2_i[2] ) - w2 * ( pl1_i[0] - pl2_i[0] ) ) / q2;
	else if ( maxq == fabs(q3) ) t = ( v2*( pl1_i[2] - pl2_i[2] ) - w2 * ( pl1_i[1] - pl2_i[1] ) ) / q3;
	else
	{
		return false;
	}

	if ( t<-( 1e-10 ) || t> dL1Length + 1e-10 ) return false;

	q1 = u2 * v1 - u1 * v2;
	q2 = u2 * w1 - u1 * w2;
	q3 = v2 * w1 - v1 * w2;

	maxq = __max(__max(fabs(q1), fabs(q2)), fabs(q3));
	if ( maxq < eps ) return false;
	if ( maxq == fabs(q1) ) t = ( u1*( pl2_i[1] - pl1_i[1] ) - v1 * ( pl2_i[0] - pl1_i[0] ) ) / q1;
	else if ( maxq == fabs(q2) ) t = ( u1*( pl2_i[2] - pl1_i[2] ) - w1 * ( pl2_i[0] - pl1_i[0] ) ) / q2;
	else if ( maxq == fabs(q3) ) t = ( v1*( pl2_i[2] - pl1_i[2] ) - w1 * ( pl2_i[1] - pl1_i[1] ) ) / q3;
	else
	{
		return false;
	}

	if ( t < -( 1e-10 ) || t > dL2Length + 1e-10 ) return false;

	// pl2선상의 교점을 구한다.
	double p2x = pl2_i[0] + u2 * t;
	double p2y = pl2_i[1] + v2 * t;
	double p2z = pl2_i[2] + w2 * t;

	pInts[0] = p2x;
	pInts[1] = p2y;
	pInts[2] = p2z;
	return true;
}

void CDgnMathFunc::mathTranUCS(const double ux[3], const double uy[3], const double uz[3], const int nData, double coor[][3])
{
	double x, y, z;
	for ( int i = 0; i < nData; i++ )
	{
		x = coor[i][0];
		y = coor[i][1];
		z = coor[i][2];
		coor[i][0] = x * ux[0] + y * ux[1] + z * ux[2];
		coor[i][1] = x * uy[0] + y * uy[1] + z * uy[2];
		coor[i][2] = x * uz[0] + y * uz[1] + z * uz[2];
	}
}

void CDgnMathFunc::mathGCS2UCS(double& x, double& y, double& z, double ucs[3][3])
{// UCS중심은 원점에 있는 것으로 가정한다. 따라서 GCS의 원하는 점은 UCS원점만큼 뺀후 입력해야한다.
	double ux = x * ucs[0][0] + y * ucs[0][1] + z * ucs[0][2];
	double uy = x * ucs[1][0] + y * ucs[1][1] + z * ucs[1][2];
	double uz = x * ucs[2][0] + y * ucs[2][1] + z * ucs[2][2];
	x = ux, y = uy, z = uz;
}

void CDgnMathFunc::mathUCS2GCS(double& x, double& y, double& z, double ucs[3][3])
{
	double gx = x * ucs[0][0] + y * ucs[1][0] + z * ucs[2][0];
	double gy = x * ucs[0][1] + y * ucs[1][1] + z * ucs[2][1];
	double gz = x * ucs[0][2] + y * ucs[1][2] + z * ucs[2][2];
	x = gx, y = gy, z = gz;
}

bool CDgnMathFunc::mathPolyCentroid(int n, double x[], double y[], double& xCentroid, double& yCentroid, double& area)
{
	xCentroid = 0, yCentroid = 0, area = 0;
	if ( n < 3 ) return false;

	int i, j;
	double ai, atmp = 0.0, xtmp = 0.0, ytmp = 0.0;
	for ( i = n - 1, j = 0; j < n; i = j, j++ )
	{
		ai = x[i] * y[j] - x[j] * y[i];
		atmp += ai;
		xtmp += ( x[j] + x[i] ) * ai;
		ytmp += ( y[j] + y[i] ) * ai;
	}
	if ( atmp == 0.0 ) return false;

	area = fabs(atmp) / 2.0;
	xCentroid = xtmp / ( 3.0 * atmp );
	yCentroid = ytmp / ( 3.0 * atmp );
	return true;
}

///// 직교 좌표계의 한 점을 원통 좌표계로 바꾸는 함수(x, y, z) => (radius, angle, height)  /////
void CDgnMathFunc::mathConvertToCylinderCoordi(double node_xyz[3], double cyn_xyz[3], double org_xyz[3], double rot_xyz[3], double pol_xyz[3])
// 대상 좌표  // 변환된 원통 좌표   // 극점좌표      // 회전축 위의 좌표  // 극축 위의 좌표
{
	double DirVec[3], nDirVec[3];  // 방향벡터
	DirVec[0] = rot_xyz[0] - org_xyz[0];
	DirVec[1] = rot_xyz[1] - org_xyz[1];
	DirVec[2] = rot_xyz[2] - org_xyz[2];
	mathNormalize(DirVec, nDirVec);
	// radius calculation.
	cyn_xyz[0] = mathDistanceFromIntersectPointToLine(org_xyz, rot_xyz, node_xyz); // 대상 좌표와 회전축과의 수직거리

	// rotating angle calculation
	double VecNode1[3], VecNode2[3], VecPol[3], nVecNode1[3], nVecNode2[3], nVecPol[3];
	double crossx, crossy, crossz;  // 대상좌표와 회전축과의 교점
	bool bLinePoint = false;  // 대상 좌표가 극축 위에 존재하는지 나타내는 flag
	// 대상좌표와 극축과의 교점구하기
	mathPLCrossPoint(org_xyz[0], org_xyz[1], org_xyz[2], nDirVec[0], nDirVec[1], nDirVec[2], node_xyz[0], node_xyz[1], node_xyz[2], crossx, crossy, crossz);
	VecNode1[0] = node_xyz[0] - crossx;
	VecNode1[1] = node_xyz[1] - crossy;
	VecNode1[2] = node_xyz[2] - crossz;
	VecPol[0] = pol_xyz[0] - org_xyz[0];
	VecPol[1] = pol_xyz[1] - org_xyz[1];
	VecPol[2] = pol_xyz[2] - org_xyz[2];
	double rotateangle = 90.0;
	mathRotate(rotateangle, nDirVec[0], nDirVec[1], nDirVec[2], node_xyz[0], node_xyz[1], node_xyz[2]);
	VecNode2[0] = node_xyz[0] - crossx;
	VecNode2[1] = node_xyz[1] - crossy;
	VecNode2[2] = node_xyz[2] - crossz;

	double angle  = 0.0;
	double angle1 = 0.0;
	double angle2 = 0.0;
	if ( !mathNormalize(VecNode1, nVecNode1) ) bLinePoint = true;
	if ( !mathNormalize(VecNode2, nVecNode2) ) bLinePoint = true;
	if ( !mathNormalize(VecPol, nVecPol) ) bLinePoint = true;
	if ( bLinePoint ) angle1 = 0;  // 대상 좌표가 극축 위에 존재할때 회전각은 0
	else
	{
		angle1 = mathCrossAngle(nVecPol, nVecNode1);  // 원래 좌표의 회전각
		angle2 = mathCrossAngle(nVecPol, nVecNode2);  // 원래 좌표를 회전축에 대해 90도 회전한 좌표의 회전각
	}
	double dTol = 0.001;
	if ( fabs(angle2 - angle1 -  90.0) < dTol ) angle = angle1;  // 1사분면
	else if ( fabs(angle2 + angle1 - 270.0) < dTol ) angle = angle1;  // 2사분면
	else if ( fabs(angle2 - angle1 +  90.0) < dTol ) angle = 360.0 - angle1; // 3사분면
	else if ( fabs(angle2 + angle1 -  90.0) < dTol ) angle = 360.0 - angle1; // 4사분면
	cyn_xyz[1] = angle;

	// height calculation
	double length = mathLength(crossx, crossy, crossz, org_xyz[0], org_xyz[1], org_xyz[2]);
	if ( fabs(crossz - org_xyz[2] - nDirVec[2] * length) < 0.0001 )
	{
		cyn_xyz[2] = length;  // 높이
	}
	else
	{
		cyn_xyz[2] = -length;
	}
}

bool CDgnMathFunc::project_on_cylinder_normal(const double node_xyz[3], const double radius, const double height, const double tolerance, double projected_node_xyz[3])
{
	if ( !( height > 0 ) ) return false;

	double norm = sqrt(pow(node_xyz[0], 2.0) + pow(node_xyz[1], 2.0));
	if ( norm < tolerance ) return false;

	double polar_angle;
	if ( fabs(node_xyz[0]) < tolerance )
	{
		polar_angle = m_pi / 2.0;
	}
	else
	{
		polar_angle = acos(node_xyz[0] / norm);
	}
	if ( node_xyz[1] < 0.0 ) polar_angle *= -1.0;

	projected_node_xyz[0] = radius * cos(polar_angle);
	projected_node_xyz[1] = radius * sin(polar_angle);
	projected_node_xyz[2] = node_xyz[2];
	return true;
}

bool CDgnMathFunc::project_on_cylinder_vector(const double node_xyz[3], const double vector_xyz[3], const double radius, const double height, const double tolerance, double projected_node_xyz[3])
{
	double line_norm = sqrt(pow(vector_xyz[0], 2.0) + pow(vector_xyz[1], 2.0) + pow(vector_xyz[2], 2.0));
	double v1 = vector_xyz[0] / line_norm;
	double v2 = vector_xyz[1] / line_norm;
	double v3 = vector_xyz[2] / line_norm;
	if ( fabs(v1) < tolerance && fabs(v2) < tolerance && fabs(v3) < 0.17452 ) return false;

	double x1 = node_xyz[0];
	double y1 = node_xyz[1];
	double z1 = node_xyz[2];
	double x2 = x1 + vector_xyz[0];
	double y2 = y1 + vector_xyz[1];
	double z2 = z1 + vector_xyz[2];
	double a11 = 1.0;
	double a22 = 1.0;
	double a33 = 0.0;
	double a34 = 0.0;
	double a44 = -( pow(radius, 2.0) );
	double A = a11 * pow(x2 - x1, 2.0) + a22 * pow(y2 - y1, 2.0) + a33 * pow(z2 - z1, 2.0);
	double B = a11 * x1 * ( x2 - x1 ) + a22 * y1 * ( y2 - y1 ) + ( a33 * z1 + a34 ) * ( z2 - z1 );
	double C = a11 * pow(x1, 2.0) + a22 * pow(y1, 2.0) + a33 * pow(z1, 2.0) + 2.0 * a34 * z1 + a44;

	if ( fabs(A) < tolerance )
	{
		if ( fabs(B) < tolerance ) return false;
		double xi = -C / ( 2.0 * B );
		projected_node_xyz[0] = x1 + xi * ( x2 - x1 );
		projected_node_xyz[1] = y1 + xi * ( y2 - y1 );
		projected_node_xyz[2] = z1 + xi * ( z2 - z1 );
		return true;
	}

	double discriminant = pow(B, 2.0) - A * C;
	if ( discriminant < 0.0 ) return false;

	double xi1 = ( -B + sqrt(discriminant) ) / A;
	double xi2 = ( -B - sqrt(discriminant) ) / A;
	double ip_xyz1[3], ip_xyz2[3];
	ip_xyz1[0] = x1 + xi1 * ( x2 - x1 );
	ip_xyz1[1] = y1 + xi1 * ( y2 - y1 );
	ip_xyz1[2] = z1 + xi1 * ( z2 - z1 );
	ip_xyz2[0] = x1 + xi2 * ( x2 - x1 );
	ip_xyz2[1] = y1 + xi2 * ( y2 - y1 );
	ip_xyz2[2] = z1 + xi2 * ( z2 - z1 );

	double distance_square1 = 0.0;
	for ( int i = 0; i < 3; i++ )
	{
		distance_square1 += pow(node_xyz[i] - ip_xyz1[i], 2.0);
	}
	double distance_square2 = 0.0;
	for ( int i = 0; i < 3; i++ )
	{
		distance_square2 += pow(node_xyz[i] - ip_xyz2[i], 2.0);
	}

	if ( distance_square1 < distance_square2 )
	{
		for ( int i = 0; i < 3; i++ )
		{
			projected_node_xyz[i] = ip_xyz1[i];
		}
		return true;
	}
	for ( int i = 0; i < 3; i++ )
	{
		projected_node_xyz[i] = ip_xyz2[i];
	}
	return true;
}

bool CDgnMathFunc::project_on_cone_normal(const double node_xyz[3], const double radius1, const double radius2, const double height, const double tolerance, double projected_node_xyz[3])
{
	if ( !( height > 0 ) ) return false;
	if ( fabs(radius1 - radius2) < tolerance ) return project_on_cylinder_normal(node_xyz, radius1, height, tolerance, projected_node_xyz);

	double norm = sqrt(pow(node_xyz[0], 2.0) + pow(node_xyz[1], 2.0));
	if ( norm < tolerance ) return false;

	double polar_angle;
	if ( fabs(node_xyz[0]) < tolerance )
	{
		polar_angle = m_pi / 2.0;
	}
	else
	{
		polar_angle = acos(node_xyz[0] / norm);
	}
	if ( node_xyz[1] < 0.0 ) polar_angle *= -1.0;

	double node_in_xz_plane_x = node_xyz[0] * cos(-polar_angle) - node_xyz[1] * sin(-polar_angle);
	double node_in_xz_plane_z = node_xyz[2];
	double line_in_xz_plane_x1 = radius1;
	double line_in_xz_plane_z1 = 0.0;
	double line_in_xz_plane_x2 = radius2;
	double line_in_xz_plane_z2 = height;
	double a = ( line_in_xz_plane_z2 - line_in_xz_plane_z1 ) / ( line_in_xz_plane_x2 - line_in_xz_plane_x1 );
	double b = line_in_xz_plane_z1 - a * line_in_xz_plane_x1;
	double c = node_in_xz_plane_z + 1.0 / a * node_in_xz_plane_x;
	double projected_node_in_xz_plane_x = ( c - b ) / ( a + 1.0 / a );
	double projected_node_in_xz_plane_z = a * projected_node_in_xz_plane_x + b;

	double xi = ( projected_node_in_xz_plane_x - line_in_xz_plane_x1 ) / ( line_in_xz_plane_x2 - line_in_xz_plane_x1 );
	if ( xi < 0.0 || xi > 1.0 ) return false;

	xi = ( projected_node_in_xz_plane_z - line_in_xz_plane_z1 ) / ( line_in_xz_plane_z2 - line_in_xz_plane_z1 );
	if ( xi < 0.0 || xi > 1.0 ) return false;

	projected_node_xyz[0] = projected_node_in_xz_plane_x * cos(polar_angle);
	projected_node_xyz[1] = projected_node_in_xz_plane_x * sin(polar_angle);
	projected_node_xyz[2] = projected_node_in_xz_plane_z;
	return true;
}

bool CDgnMathFunc::project_on_cone_vector(const double node_xyz[3], const double vector_xyz[3], const double radius1, const double radius2, const double height, const double tolerance, double projected_node_xyz[3])
{
	double x1 = node_xyz[0];
	double y1 = node_xyz[1];
	double z1 = node_xyz[2];
	double x2 = x1 + vector_xyz[0];
	double y2 = y1 + vector_xyz[1];
	double z2 = z1 + vector_xyz[2];
	double a11 = 1.0;
	double a22 = 1.0;
	double a33 = -pow(( radius2 - radius1 ) / height, 2.0);
	double a34 = -2.0 * radius1 * ( radius2 - radius1 ) / height;
	double a44 = -pow(radius1, 2.0);
	double A = a11 * pow(x2 - x1, 2.0) + a22 * pow(y2 - y1, 2.0) + a33 * pow(z2 - z1, 2.0);
	double B = a11 * x1 * ( x2 - x1 ) + a22 * y1 * ( y2 - y1 ) + ( a33 * z1 + a34 ) * ( z2 - z1 );
	double C = a11 * pow(x1, 2.0) + a22 * pow(y1, 2.0) + a33 * pow(z1, 2.0) + 2.0 * a34 * z1 + a44;

	if ( fabs(A) < tolerance )
	{
		if ( fabs(B) < tolerance ) return false;
		double xi = -C / ( 2.0 * B );
		projected_node_xyz[0] = x1 + xi * ( x2 - x1 );
		projected_node_xyz[1] = y1 + xi * ( y2 - y1 );
		projected_node_xyz[2] = z1 + xi * ( z2 - z1 );
		return true;
	}

	double discriminant = pow(B, 2.0) - A * C;
	if ( discriminant < 0.0 ) return false;
	double xi1 = ( -B + sqrt(discriminant) ) / A;
	double xi2 = ( -B - sqrt(discriminant) ) / A;
	double ip_xyz1[3], ip_xyz2[3];
	ip_xyz1[0] = x1 + xi1 * ( x2 - x1 );
	ip_xyz1[1] = y1 + xi1 * ( y2 - y1 );
	ip_xyz1[2] = z1 + xi1 * ( z2 - z1 );
	ip_xyz2[0] = x1 + xi2 * ( x2 - x1 );
	ip_xyz2[1] = y1 + xi2 * ( y2 - y1 );
	ip_xyz2[2] = z1 + xi2 * ( z2 - z1 );

	double distance_square1 = 0.0;
	for ( int i = 0; i < 3; i++ )
	{
		distance_square1 += pow(node_xyz[i] - ip_xyz1[i], 2.0);
	}
	double distance_square2 = 0.0;
	for ( int i = 0; i < 3; i++ )
	{
		distance_square2 += pow(node_xyz[i] - ip_xyz2[i], 2.0);
	}

	if ( distance_square1 < distance_square2 )
	{
		for ( int i = 0; i < 3; i++ )
		{
			projected_node_xyz[i] = ip_xyz1[i];
		}
		return true;
	}
	for ( int i = 0; i < 3; i++ )
	{
		projected_node_xyz[i] = ip_xyz2[i];
	}
	return true;
}

bool CDgnMathFunc::project_on_sphere_normal(const double node_xyz[3], const double radius, const double tolerance, double projected_node_xyz[3])
{
	if ( !( radius > 0 ) ) return false;

	double norm = sqrt(pow(node_xyz[0], 2.0) + pow(node_xyz[1], 2.0));
	if ( norm < tolerance )
	{
		if ( fabs(node_xyz[2]) < tolerance ) return false;
		projected_node_xyz[0] = projected_node_xyz[1] = 0.0;
		if ( node_xyz[2] > 0.0 ) projected_node_xyz[2] = radius;
		else projected_node_xyz[2] = -radius;
		return true;
	}

	double polar_angle;
	if ( fabs(node_xyz[0]) < tolerance )
	{
		polar_angle = m_pi / 2.0;
	}
	else
	{
		polar_angle = acos(node_xyz[0] / norm);
	}
	if ( node_xyz[1] < 0.0 ) polar_angle *= -1.0;

	double node_in_xz_plane_x = node_xyz[0] * cos(-polar_angle) - node_xyz[1] * sin(-polar_angle);
	double node_in_xz_plane_z = node_xyz[2];
	if ( fabs(node_in_xz_plane_z) < 0.0 )
	{
		projected_node_xyz[0] = radius * cos(polar_angle);
		projected_node_xyz[1] = radius * sin(polar_angle);
		projected_node_xyz[2] = 0.0;
		return true;
	}
	double projected_node_in_xz_plane_x, projected_node_in_xz_plane_z;
	if ( node_in_xz_plane_z > 0.0 )
	{
		projected_node_in_xz_plane_z = sqrt(pow(radius, 2.0) / ( 1.0 + pow(node_in_xz_plane_x / node_in_xz_plane_z, 2.0) ));
	}
	else
	{
		projected_node_in_xz_plane_z = -sqrt(pow(radius, 2.0) / ( 1.0 + pow(node_in_xz_plane_x / node_in_xz_plane_z, 2.0) ));
	}
	projected_node_in_xz_plane_x = sqrt(pow(radius, 2.0) - pow(projected_node_in_xz_plane_z, 2.0));

	projected_node_xyz[0] = projected_node_in_xz_plane_x * cos(polar_angle);
	projected_node_xyz[1] = projected_node_in_xz_plane_x * sin(polar_angle);
	projected_node_xyz[2] = projected_node_in_xz_plane_z;
	return true;
}

bool CDgnMathFunc::project_on_sphere_vector(const double node_xyz[3], const double vector_xyz[3], const double radius, const double tolerance, double projected_node_xyz[3])
{
	double x1 = node_xyz[0];
	double y1 = node_xyz[1];
	double z1 = node_xyz[2];
	double x2 = x1 + vector_xyz[0];
	double y2 = y1 + vector_xyz[1];
	double z2 = z1 + vector_xyz[2];
	double a11 = 1.0;
	double a22 = 1.0;
	double a33 = 1.0;
	double a34 = 0.0;
	double a44 = -( pow(radius, 2.0) );
	double A = a11 * pow(x2 - x1, 2.0) + a22 * pow(y2 - y1, 2.0) + a33 * pow(z2 - z1, 2.0);
	double B = a11 * x1 * ( x2 - x1 ) + a22 * y1 * ( y2 - y1 ) + ( a33 * z1 + a34 ) * ( z2 - z1 );
	double C = a11 * pow(x1, 2.0) + a22 * pow(y1, 2.0) + a33 * pow(z1, 2.0) + 2.0 * a34 * z1 + a44;

	if ( fabs(A) < tolerance )
	{
		if ( fabs(B) < tolerance ) return false;
		double xi = -C / ( 2.0 * B );
		projected_node_xyz[0] = x1 + xi * ( x2 - x1 );
		projected_node_xyz[1] = y1 + xi * ( y2 - y1 );
		projected_node_xyz[2] = z1 + xi * ( z2 - z1 );
		return true;
	}

	double discriminant = pow(B, 2.0) - A * C;
	if ( discriminant < 0.0 ) return false;

	double xi1 = ( -B + sqrt(discriminant) ) / A;
	double xi2 = ( -B - sqrt(discriminant) ) / A;
	double ip_xyz1[3], ip_xyz2[3];
	ip_xyz1[0] = x1 + xi1 * ( x2 - x1 );
	ip_xyz1[1] = y1 + xi1 * ( y2 - y1 );
	ip_xyz1[2] = z1 + xi1 * ( z2 - z1 );
	ip_xyz2[0] = x1 + xi2 * ( x2 - x1 );
	ip_xyz2[1] = y1 + xi2 * ( y2 - y1 );
	ip_xyz2[2] = z1 + xi2 * ( z2 - z1 );

	double distance_square1 = 0.0;
	for ( int i = 0; i < 3; i++ )
	{
		distance_square1 += pow(node_xyz[i] - ip_xyz1[i], 2.0);
	}
	double distance_square2 = 0.0;
	for ( int i = 0; i < 3; i++ )
	{
		distance_square2 += pow(node_xyz[i] - ip_xyz2[i], 2.0);
	}

	if ( distance_square1 < distance_square2 )
	{
		for ( int i = 0; i < 3; i++ )
		{
			projected_node_xyz[i] = ip_xyz1[i];
		}
		return true;
	}
	for ( int i = 0; i < 3; i++ )
	{
		projected_node_xyz[i] = ip_xyz2[i];
	}
	return true;
}

bool CDgnMathFunc::project_on_ellipsoid_normal(const double node_xyz[3], const double radius, const double height, const double tolerance, double projected_node_xyz[3])
{
	if ( !( radius > 0 && height > 0.0 ) ) return false;

	double norm = sqrt(pow(node_xyz[0], 2.0) + pow(node_xyz[1], 2.0));
	if ( norm < tolerance )
	{
		projected_node_xyz[0] = projected_node_xyz[1] = 0.0;
		projected_node_xyz[2] = height;
		return true;
	}

	double polar_angle;
	if ( fabs(node_xyz[0]) < tolerance )
	{
		polar_angle = m_pi / 2.0;
	}
	else
	{
		polar_angle = acos(node_xyz[0] / norm);
	}
	if ( node_xyz[1] < 0.0 ) polar_angle *= -1.0;

	double node_in_xz_plane_x = node_xyz[0] * cos(-polar_angle) - node_xyz[1] * sin(-polar_angle);
	double node_in_xz_plane_z = node_xyz[2];
	if ( fabs(node_in_xz_plane_z) < 0.0 )
	{
		if ( node_in_xz_plane_x > radius ) return false;
		projected_node_xyz[0] = radius * cos(polar_angle);
		projected_node_xyz[1] = radius * sin(polar_angle);
		projected_node_xyz[2] = 0.0;
		return true;
	}

	double projected_node_in_xz_plane_x, projected_node_in_xz_plane_z, m;
	if ( radius > height )
	{
		double focus1 = sqrt(pow(radius, 2.0) - pow(height, 2.0));
		double focus2 = -focus1;
		double norm1 = sqrt(pow(node_in_xz_plane_x - focus1, 2.0) + pow(node_in_xz_plane_z, 2.0));
		double norm2 = sqrt(pow(node_in_xz_plane_x - focus2, 2.0) + pow(node_in_xz_plane_z, 2.0));
		double angle1 = acos(( node_in_xz_plane_x - focus1 ) / norm1);
		double angle2 = acos(( node_in_xz_plane_x - focus2 ) / norm2);
		m = tan(( angle1 + angle2 ) / 2.0);
	}
	else
	{
		double focus1 = sqrt(pow(height, 2.0) - pow(radius, 2.0));
		double focus2 = -focus1;
		double norm1 = sqrt(pow(node_in_xz_plane_x, 2.0) + pow(node_in_xz_plane_z - focus1, 2.0));
		double norm2 = sqrt(pow(node_in_xz_plane_x, 2.0) + pow(node_in_xz_plane_z - focus2, 2.0));
		double angle1 = acos(( node_in_xz_plane_z - focus1 ) / norm1);
		double angle2 = acos(( node_in_xz_plane_z - focus2 ) / norm2);
		m = tan(m_pi / 2.0 - ( angle1 + angle2 ) / 2.0);
	}
	projected_node_in_xz_plane_x = radius / sqrt(1.0 + pow(height / radius * m, 2.0));
	projected_node_in_xz_plane_z = pow(height / radius, 2.0) * m * projected_node_in_xz_plane_x;

	projected_node_xyz[0] = projected_node_in_xz_plane_x * cos(polar_angle);
	projected_node_xyz[1] = projected_node_in_xz_plane_x * sin(polar_angle);
	projected_node_xyz[2] = projected_node_in_xz_plane_z;
	return true;
}

bool CDgnMathFunc::project_on_ellipsoid_vector(const double node_xyz[3], const double vector_xyz[3], const double radius, const double height, const double tolerance, double projected_node_xyz[3])
{
	double x1 = node_xyz[0];
	double y1 = node_xyz[1];
	double z1 = node_xyz[2];
	double x2 = x1 + vector_xyz[0];
	double y2 = y1 + vector_xyz[1];
	double z2 = z1 + vector_xyz[2];
	double a11 = 1.0 / pow(radius, 2.0);
	double a22 = 1.0 / pow(radius, 2.0);
	double a33 = 1.0 / pow(height, 2.0);
	double a34 = 0.0;
	double a44 = -1.0;
	double A = a11 * pow(x2 - x1, 2.0) + a22 * pow(y2 - y1, 2.0) + a33 * pow(z2 - z1, 2.0);
	double B = a11 * x1 * ( x2 - x1 ) + a22 * y1 * ( y2 - y1 ) + ( a33 * z1 + a34 ) * ( z2 - z1 );
	double C = a11 * pow(x1, 2.0) + a22 * pow(y1, 2.0) + a33 * pow(z1, 2.0) + 2.0 * a34 * z1 + a44;

	if ( fabs(A) < tolerance )
	{
		if ( fabs(B) < tolerance ) return false;
		double xi = -C / ( 2.0 * B );
		projected_node_xyz[0] = x1 + xi * ( x2 - x1 );
		projected_node_xyz[1] = y1 + xi * ( y2 - y1 );
		projected_node_xyz[2] = z1 + xi * ( z2 - z1 );
		return true;
	}

	double discriminant = pow(B, 2.0) - A * C;
	if ( discriminant < 0.0 ) return false;
	double xi1 = ( -B + sqrt(discriminant) ) / A;
	double xi2 = ( -B - sqrt(discriminant) ) / A;
	double ip_xyz1[3], ip_xyz2[3];
	ip_xyz1[0] = x1 + xi1 * ( x2 - x1 );
	ip_xyz1[1] = y1 + xi1 * ( y2 - y1 );
	ip_xyz1[2] = z1 + xi1 * ( z2 - z1 );
	ip_xyz2[0] = x1 + xi2 * ( x2 - x1 );
	ip_xyz2[1] = y1 + xi2 * ( y2 - y1 );
	ip_xyz2[2] = z1 + xi2 * ( z2 - z1 );

	double distance_square1 = 0.0;
	for ( int i = 0; i < 3; i++ )
	{
		distance_square1 += pow(node_xyz[i] - ip_xyz1[i], 2.0);
	}
	double distance_square2 = 0.0;
	for ( int i = 0; i < 3; i++ )
	{
		distance_square2 += pow(node_xyz[i] - ip_xyz2[i], 2.0);
	}

	if ( distance_square1 < distance_square2 )
	{
		for ( int i = 0; i < 3; i++ )
		{
			projected_node_xyz[i] = ip_xyz1[i];
		}
		return true;
	}
	for ( int i = 0; i < 3; i++ )
	{
		projected_node_xyz[i] = ip_xyz2[i];
	}
	return true;
}

bool CDgnMathFunc::project_on_plane_normal(const double node_xyz[3], const double plane_xyz1[3], const double plane_xyz2[3], const double plane_xyz3[3], const double tolerance, double projected_node_xyz[3])
{
	double vector12[3], vector13[3];
	for ( int i = 0; i < 3; i++ )
	{
		vector12[i] = plane_xyz2[i] - plane_xyz1[i];
		vector13[i] = plane_xyz3[i] - plane_xyz1[i];
	}
	double a = vector12[1] * vector13[2] - vector13[1] * vector12[2];
	double b = vector13[0] * vector12[2] - vector12[0] * vector13[2];
	double c = vector12[0] * vector13[1] - vector13[0] * vector12[1];
	double d = a * plane_xyz1[0] + b * plane_xyz1[1] + c * plane_xyz1[2];
	double denominator = pow(a, 2.0) + pow(b, 2.0) + pow(c, 2.0);
	if ( fabs(denominator) < tolerance ) return false;

	double xi = ( d - a * node_xyz[0] - b * node_xyz[1] - c * node_xyz[2] ) / denominator;
	projected_node_xyz[0] = node_xyz[0] + xi * a;
	projected_node_xyz[1] = node_xyz[1] + xi * b;
	projected_node_xyz[2] = node_xyz[2] + xi * c;
	return true;
}

bool CDgnMathFunc::project_on_plane_vector(const double node_xyz[3], const double vector_xyz[3], const double plane_xyz1[3], const double plane_xyz2[3], const double plane_xyz3[3], double projected_node_xyz[3])
{
	return project_on_plane_vector(node_xyz, vector_xyz, plane_xyz1, plane_xyz2, plane_xyz3, projected_node_xyz, 0.0174524);
}

bool CDgnMathFunc::project_on_plane_vector_for_CuttingPlane(const double node_xyz[3], const double vector_xyz[3], const double plane_xyz1[3], const double plane_xyz2[3], const double plane_xyz3[3], double projected_node_xyz[3])
{
	return project_on_plane_vector(node_xyz, vector_xyz, plane_xyz1, plane_xyz2, plane_xyz3, projected_node_xyz, 1e-10);
}

bool CDgnMathFunc::project_on_plane_vector(const double node_xyz[3], const double vector_xyz[3], const double plane_xyz1[3], const double plane_xyz2[3], const double plane_xyz3[3], double projected_node_xyz[3], double tolerance)
{
	double vector12[3], vector13[3];
	for ( int i = 0; i < 3; i++ )
	{
		vector12[i] = plane_xyz2[i] - plane_xyz1[i];
		vector13[i] = plane_xyz3[i] - plane_xyz1[i];
	}
	double a = vector12[1] * vector13[2] - vector13[1] * vector12[2];
	double b = vector13[0] * vector12[2] - vector12[0] * vector13[2];
	double c = vector12[0] * vector13[1] - vector13[0] * vector12[1];
	double plane_norm = sqrt(pow(a, 2.0) + pow(b, 2.0) + pow(c, 2.0));
	if ( (plane_norm == 0.0) ) return false; // not plane

	a /= plane_norm;
	b /= plane_norm;
	c /= plane_norm;
	double d = a * plane_xyz1[0] + b * plane_xyz1[1] + c * plane_xyz1[2];
	double line_norm = sqrt(pow(vector_xyz[0], 2.0) + pow(vector_xyz[1], 2.0) + pow(vector_xyz[2], 2.0));
	double v1 = vector_xyz[0] / line_norm;
	double v2 = vector_xyz[1] / line_norm;
	double v3 = vector_xyz[2] / line_norm;
	double denominator = a * v1 + b * v2 + c * v3;
	if ( fabs(denominator) < tolerance ) return false;

	double xi = ( d - a * node_xyz[0] - b * node_xyz[1] - c * node_xyz[2] ) / denominator;
	projected_node_xyz[0] = node_xyz[0] + xi * v1;
	projected_node_xyz[1] = node_xyz[1] + xi * v2;
	projected_node_xyz[2] = node_xyz[2] + xi * v3;
	return true;
}

bool CDgnMathFunc::check_quad_curve(const double curve_xyz_1[3], const double curve_xyz_2[3], const double curve_xyz_3[3], const double ZERO)
{
	double matrix[3][3];
	matrix[0][0] = pow(curve_xyz_1[0], 2.0);
	matrix[1][0] = pow(curve_xyz_2[0], 2.0);
	matrix[2][0] = pow(curve_xyz_3[0], 2.0);
	matrix[0][1] = curve_xyz_1[0];
	matrix[1][1] = curve_xyz_2[0];
	matrix[2][1] = curve_xyz_3[0];
	matrix[0][2] = 1.0;
	matrix[1][2] = 1.0;
	matrix[2][2] = 1.0;

	double determinant =
		- matrix[0][2] * matrix[1][1] * matrix[2][0] + matrix[0][1] * matrix[1][2] * matrix[2][0]
		+ matrix[0][2] * matrix[1][0] * matrix[2][1] - matrix[0][0] * matrix[1][2] * matrix[2][1]
		- matrix[0][1] * matrix[1][0] * matrix[2][2] + matrix[0][0] * matrix[1][1] * matrix[2][2];

	if ( fabs(determinant) < ZERO ) return 0; // cannot define the quadratic curve
	return 1;
}

bool CDgnMathFunc::project_on_quad_curve_normal(const double curve_xyz_1[3], const double curve_xyz_2[3], const double curve_xyz_3[3], const double point[3], const double ZERO, double project_point[3])
// curve_xyz_1, curve_xyz_2, curve_xyz_3: 3 coordinates defining the quadratic curve (in UCS)
// point                                : coordinate of the target point (in UCS)
// ZERO                                 : zero tolerance
// project_point                        : projected coordinate (in UCS)
{
	double matrix[3][3], inverse_matrix[3][3];
	matrix[0][0] = pow(curve_xyz_1[0], 2.0);
	matrix[1][0] = pow(curve_xyz_2[0], 2.0);
	matrix[2][0] = pow(curve_xyz_3[0], 2.0);
	matrix[0][1] = curve_xyz_1[0];
	matrix[1][1] = curve_xyz_2[0];
	matrix[2][1] = curve_xyz_3[0];
	matrix[0][2] = 1.0;
	matrix[1][2] = 1.0;
	matrix[2][2] = 1.0;

	double determinant =
		-matrix[0][2] * matrix[1][1] * matrix[2][0] + matrix[0][1] * matrix[1][2] * matrix[2][0]
		+ matrix[0][2] * matrix[1][0] * matrix[2][1] - matrix[0][0] * matrix[1][2] * matrix[2][1]
		- matrix[0][1] * matrix[1][0] * matrix[2][2] + matrix[0][0] * matrix[1][1] * matrix[2][2];

	if ( fabs(determinant) < ZERO ) return 0; // cannot define the quadratic curve

	inverse_matrix[0][0] = ( -matrix[1][2] * matrix[2][1] + matrix[1][1] * matrix[2][2] ) / determinant;
	inverse_matrix[0][1] = (  matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2] ) / determinant;
	inverse_matrix[0][2] = ( -matrix[0][2] * matrix[1][1] + matrix[0][1] * matrix[1][2] ) / determinant;
	inverse_matrix[1][0] = (  matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2] ) / determinant;
	inverse_matrix[1][1] = ( -matrix[0][2] * matrix[2][0] + matrix[0][0] * matrix[2][2] ) / determinant;
	inverse_matrix[1][2] = (  matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2] ) / determinant;
	inverse_matrix[2][0] = ( -matrix[1][1] * matrix[2][0] + matrix[1][0] * matrix[2][1] ) / determinant;
	inverse_matrix[2][1] = (  matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1] ) / determinant;
	inverse_matrix[2][2] = ( -matrix[0][1] * matrix[1][0] + matrix[0][0] * matrix[1][1] ) / determinant;

	double vector_y[3], vector_abc[3];
	vector_y[0] = curve_xyz_1[1];
	vector_y[1] = curve_xyz_2[1];
	vector_y[2] = curve_xyz_3[1];
	for ( int i = 0; i < 3; i++ )
	{
		vector_abc[i] = 0.0;
		for ( int j = 0; j < 3; j++ )
		{
			vector_abc[i] += inverse_matrix[i][j] * vector_y[j];
		}
	}
	// quadratic curve: y = a x^2 + b x + c
	double a = vector_abc[0];
	double b = vector_abc[1];
	double c = vector_abc[2];

	double                        max_x = curve_xyz_1[0];
	if ( curve_xyz_2[0] > max_x ) max_x = curve_xyz_2[0];
	if ( curve_xyz_3[0] > max_x ) max_x = curve_xyz_3[0];
	double                        min_x = curve_xyz_1[0];
	if ( curve_xyz_2[0] < min_x ) min_x = curve_xyz_2[0];
	if ( curve_xyz_3[0] < min_x ) min_x = curve_xyz_3[0];
	int N_division = 100;
	double dx = ( max_x - min_x ) / N_division;
	double minimum = 1.0e+35;
	for ( int i = 0; i <= N_division; i++ )
	{
		double x = min_x + dx * i;
		double dy = 2.0 * a * x + b;
		if ( fabs(dy) < ZERO )
		{
			double y = a * pow(x, 2.0) + b * x + c;
			double distance = pow(point[0] - x, 2.0) + pow(point[1] - y, 2.0);
			if ( distance < minimum )
			{
				minimum = distance;
				project_point[0] = x;
				project_point[1] = y;
			}
			continue;
		}
		// line: y = p x + q
		double p = -1.0 / dy;
		double q = point[1] - p * point[0];
		double r = a;
		double s = b - p;
		double t = c - q;
		double D = pow(s, 2.0) - 4.0 * r * t;
		if ( D < ZERO ) continue; // no real solution

		double x1 = ( -s + sqrt(D) ) / ( 2.0 * r );
		double x2 = ( -s - sqrt(D) ) / ( 2.0 * r );
		double y1 = p * x1 + q;
		double y2 = p * x2 + q;
		double distance1 = pow(point[0] - x1, 2.0) + pow(point[1] - y1, 2.0);
		double distance2 = pow(point[0] - x2, 2.0) + pow(point[1] - y2, 2.0);
		if ( distance1 < distance2 )
		{
			if ( distance1 < minimum )
			{
				minimum = distance1;
				project_point[0] = x1;
				project_point[1] = y1;
			}
		}
		else
		{
			if ( distance2 < minimum )
			{
				minimum = distance2;
				project_point[0] = x2;
				project_point[1] = y2;
			}
		}
	}
	project_point[2] = point[2];

	return 1;
}

bool CDgnMathFunc::project_on_quad_curve_vector(const double curve_xyz_1[3], const double curve_xyz_2[3], const double curve_xyz_3[3], const double point[3], const double direction_vector[3], const double ZERO, double project_point[3])
// curve_xyz_1, curve_xyz_2, curve_xyz_3: 3 coordinates defining the quadratic curve (in UCS)
// point                                : coordinate of the target point (in UCS)
// direction_vector                     : direction vector
// ZERO                                 : zero tolerance
// project_point                        : projected coordinate (in UCS)
{
	double magnitude = pow(direction_vector[0], 2.0) + pow(direction_vector[1], 2.0) + pow(direction_vector[2], 2.0);
	if ( fabs(magnitude) < ZERO ) return 0;

	double matrix[3][3], inverse_matrix[3][3];
	matrix[0][0] = pow(curve_xyz_1[0], 2.0);
	matrix[1][0] = pow(curve_xyz_2[0], 2.0);
	matrix[2][0] = pow(curve_xyz_3[0], 2.0);
	matrix[0][1] = curve_xyz_1[0];
	matrix[1][1] = curve_xyz_2[0];
	matrix[2][1] = curve_xyz_3[0];
	matrix[0][2] = 1.0;
	matrix[1][2] = 1.0;
	matrix[2][2] = 1.0;

	double determinant =
		- matrix[0][2] * matrix[1][1] * matrix[2][0] + matrix[0][1] * matrix[1][2] * matrix[2][0]
		+ matrix[0][2] * matrix[1][0] * matrix[2][1] - matrix[0][0] * matrix[1][2] * matrix[2][1]
		- matrix[0][1] * matrix[1][0] * matrix[2][2] + matrix[0][0] * matrix[1][1] * matrix[2][2];

	if ( fabs(determinant) < ZERO ) return 0; // cannot define the quadratic curve

	inverse_matrix[0][0] = ( -matrix[1][2] * matrix[2][1] + matrix[1][1] * matrix[2][2] ) / determinant;
	inverse_matrix[0][1] = (  matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2] ) / determinant;
	inverse_matrix[0][2] = ( -matrix[0][2] * matrix[1][1] + matrix[0][1] * matrix[1][2] ) / determinant;
	inverse_matrix[1][0] = (  matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2] ) / determinant;
	inverse_matrix[1][1] = ( -matrix[0][2] * matrix[2][0] + matrix[0][0] * matrix[2][2] ) / determinant;
	inverse_matrix[1][2] = (  matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2] ) / determinant;
	inverse_matrix[2][0] = ( -matrix[1][1] * matrix[2][0] + matrix[1][0] * matrix[2][1] ) / determinant;
	inverse_matrix[2][1] = (  matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1] ) / determinant;
	inverse_matrix[2][2] = ( -matrix[0][1] * matrix[1][0] + matrix[0][0] * matrix[1][1] ) / determinant;

	double vector_y[3], vector_abc[3];
	vector_y[0] = curve_xyz_1[1];
	vector_y[1] = curve_xyz_2[1];
	vector_y[2] = curve_xyz_3[1];
	for ( int i = 0; i < 3; i++ )
	{
		vector_abc[i] = 0.0;
		for ( int j = 0; j < 3; j++ )
		{
			vector_abc[i] += inverse_matrix[i][j] * vector_y[j];
		}
	}
	// quadratic curve: y = a x^2 + b x + c
	double a = vector_abc[0];
	double b = vector_abc[1];
	double c = vector_abc[2];

	if ( fabs(direction_vector[0]) < ZERO )
	{
		project_point[0] = point[0];
		project_point[1] = a * pow(point[0], 2.0) + b * point[0] + c;
		project_point[2] = point[2];
		return 1;
	}

	// line: y = p x + q
	double p = direction_vector[1] / direction_vector[0];
	double q = point[1] - p * point[0];
	double r = a;
	double s = b - p;
	double t = c - q;
	double D = pow(s, 2.0) - 4.0 * r * t;
	if ( D < ZERO ) return 0; // no real solution

	double x1 = ( -s + sqrt(D) ) / ( 2.0 * r );
	double x2 = ( -s - sqrt(D) ) / ( 2.0 * r );
	double y1 = p * x1 + q;
	double y2 = p * x2 + q;
	double distance1 = pow(point[0] - x1, 2.0) + pow(point[1] - y1, 2.0);
	double distance2 = pow(point[0] - x2, 2.0) + pow(point[1] - y2, 2.0);
	if ( distance1 < distance2 )
	{
		project_point[0] = x1;
		project_point[1] = y1;
	}
	else
	{
		project_point[0] = x2;
		project_point[1] = y2;
	}
	project_point[2] = point[2];

	return 1;
}

double CDgnMathFunc::InterpolateLinear(const std::vector<std::pair<double, double>>& vdXY, const double& dX)
{
	const size_t nSize = vdXY.size();
	if ( nSize == 0 ) { return 0.0; }
	if ( nSize == 1 ) { return vdXY[0].second; }

	if ( dX < vdXY[0].first ) { return vdXY[0].second; }
	if ( dX > vdXY[nSize-1].first ) { return vdXY[nSize-1].second; }

	for ( int i = 0; i < nSize - 1; ++i )
	{
		const double dX1 = vdXY[i].first;
		const double dY1 = vdXY[i].second;

		const double dX2 = vdXY[i+1].first;
		const double dY2 = vdXY[i+1].second;

		if ( ( dX > dX1 || (dX == dX1) ) && ( dX < dX2 || (dX == dX2) ) )
		{
			return InterpolateLinear(dX1, dX2, dY1, dY2, dX);
		}
	}
	return vdXY[0].second;
}

double CDgnMathFunc::InterpolateLinear(const double& dX1, const double& dX2, const double& dY1, const double& dY2, const double& dX)
{
	// if ( ( dX < dX1 ) || EQ(dX, dX1) ) { return dY1; }
	// if ( ( dX > dX2 ) || EQ(dX, dX2) ) { return dY2; }
	// return dY1 + dgn::SafeDiv(dY2 - dY1, dX2 - dX1) * ( dX - dX1 );
    return 0.0;
}

///////////////////////////////////////////////////////////////////////////////
// 기능 : 두 점의 값으로부터 원하는 위치에서의 직선보간값을 구한다.(1차원 직선보간)
// 입력 : 두 점의 값 및 왼쪽점의 값으로부터 오른쪽으로의 비율(0~1)
//    dLeft                         :       dRigth
//       |<-- dDistRatioFromLeft -->|
//
double CDgnMathFunc::mathInterpolate(double dLeft, double dRight, double dDistRatioFromLeft)
{
	double dValue;
	if ( IsNormalZero(dDistRatioFromLeft) )
	{
		dValue = dLeft;
	}
	else if ( dDistRatioFromLeft > 1.0 )
	{
		dValue = dRight;
	}
	else
	{
		dValue = dLeft + ( dRight - dLeft ) * dDistRatioFromLeft;
	}

	return dValue;
}

///////////////////////////////////////////////////////////////////////////////
// 기능 : 네 모퉁이의 값으로부터 원하는 위치에서의 직선보간값을 구한다.(2차원 직선보간)
// 입력 : 네 모퉁이의 값 및 왼쪽 아래 값으로부터 오른쪽과 위쪽으로의 비율(0~1)
//
//    dLeftTop                              dRightTop
//                                                         ___
//        --------------------------*                       |
//                                  :                       |
//                                  :                 dDistRatioFromBot
//                                  :                       |
//    dLeftBot                      :       dRigthTop      ---
//       |<-- dDistRatioFromLeft -->|
//
double CDgnMathFunc::mathInterpolate(double dLeftBot, double dRightBot, double dLeftTop, double dRightTop,
	double dDistRatioFromLeft, double dDistRatioFromBot)
{
	double dValueBot, dValueTop, dValue;
	if ( IsNormalZero(dDistRatioFromLeft) )
	{
		dValueBot = dLeftBot;
		dValueTop = dLeftTop;
	}
	else if ( dDistRatioFromLeft > 1.0 )
	{
		dValueBot = dRightBot;
		dValueTop = dRightTop;
	}
	else
	{
		dValueBot = dLeftBot + ( dRightBot - dLeftBot ) * dDistRatioFromLeft;
		dValueTop = dLeftTop + ( dRightTop - dLeftTop ) * dDistRatioFromLeft;
	}

	if ( IsNormalZero(dDistRatioFromBot) )
	{
		dValue = dValueBot;
	}
	else if ( dDistRatioFromBot > 1.0 )
	{
		dValue = dValueTop;
	}
	else
	{
		dValue = dValueBot + ( dValueTop - dValueBot ) * dDistRatioFromBot;
	}

	return dValue;
}

void CDgnMathFunc::mathCrossPoint2D(double line_i[2], double line_j[2], double point[2], double point_cross[2])
{
	double a, b, l, m, x, y;
	double dx, dy;

	x = point[0], y = point[1];
	a = line_i[0], b = line_i[1];
	dx = line_j[0] - a, dy = line_j[1] - b;

	double dLength = mathLength(dx, dy);
	mathNormalize(dx, dy, l, m);
	double t;
	if ( IsNormalZero(l) && IsNormalZero(m) )
	{
		t = 0.0;
	}
	else
	{
		t = ( ( x - a )*l + ( y - b )*m ) / ( l*l + m * m );
	}

	point_cross[0] = a + t * l;
	point_cross[1] = b + t * m;
}

double CDgnMathFunc::mathDistanceFromIntersectPointToLine2D(double line_i[2], double line_j[2], double point[2], int& nState)
{
	double a, b, l, m, x, y;
	double dx, dy, px, py;

	x = point[0], y = point[1];
	a = line_i[0], b = line_i[1];
	dx = line_j[0] - a, dy = line_j[1] - b;

	double dLength = mathLength(dx, dy);
	mathNormalize(dx, dy, l, m);
	double t;
	if ( IsNormalZero(l) && IsNormalZero(m) )
	{
		t = 0.0;
	}
	else
	{
		t = ( ( x - a )*l + ( y - b )*m ) / ( l*l + m * m );
	}

	if ( t < 0. ) { nState = 0; }
	else if ( t > dLength ) { nState = 2; }
	else { nState = 1; }

	px = a + t * l;
	py = b + t * m;

	return mathLength(point[0] - px, point[1] - py);
}

double CDgnMathFunc::mathCrossAngle2DSign(double u1x, double u1y, double u2x, double u2y)
{// Normalized 된 2차원 벡터(z=0)만 써야한다.
	double vector13D[3], vector23D[3];
	vector13D[0] = u1x;  vector13D[1] = u1y; vector13D[2] = 0.0;
	vector23D[0] = u2x;  vector23D[1] = u2y; vector23D[2] = 0.0;

	return mathCrossAngle2D(vector13D, vector23D);
}

double CDgnMathFunc::mathCross2D(double vector1[2], double vector2[2])
{
	return vector1[0] * vector2[1] - vector1[1] * vector2[0];
}

bool CDgnMathFunc::mathNormalize2D(double vector[2], double vectorn[2])
{
	return mathNormalize(vector[0], vector[1], vectorn[0], vectorn[1]);
}

void CDgnMathFunc::mathSwap2D(double p1[2], double p2[2])
{
	double tx = p1[0];
	double ty = p1[1];

	p1[0] = p2[0];
	p1[1] = p2[1];

	p2[0] = tx;
	p2[1] = ty;
}

bool CDgnMathFunc::mathIsPointOfLine2D(double bound1[2], double bound2[2], double targetPt[2], bool isOnLine = true)
{// SHIN (2006.2.4) e-15로 했을경우 오차범위가 너무 좁음  // SHIN (2006.2.28) 오차범위를 선분의 길이의 e-15배로 수정// SHIN (2006.6.28) 오차범위를 선분의 길이의 e-11배로 수정
	double dNormalZero = m_NormalZero * 10000.0;
	double dTol = __max(mathLength(bound1[0] - bound2[0], bound1[1] - bound2[1]), m_NormalZero);
	dTol = __min(dTol, 1.0 / dNormalZero); // MQC.10619 3D PM-Curve에서 강도점을 찾지 못함. dNormalZero*dTol < 1.0 이 되도록 수정.

	if ( fabs(targetPt[0] - bound1[0]) < dNormalZero*dTol && fabs(targetPt[1] - bound1[1]) < dNormalZero*dTol )
	{
		return isOnLine;
	}
	if ( fabs(targetPt[0] - bound2[0]) < dNormalZero*dTol && fabs(targetPt[1] - bound2[1]) < dNormalZero*dTol )
	{
		return isOnLine;
	}

	double vector1[2], vector2[2];
	vector1[0] = bound1[0] - targetPt[0];
	vector1[1] = bound1[1] - targetPt[1];
	mathNormalize2D(vector1, vector1);
	vector2[0] = bound2[0] - targetPt[0];
	vector2[1] = bound2[1] - targetPt[1];
	mathNormalize2D(vector2, vector2);

	// 외적(vector product , ^ )이 0이면 두 벡터는 평행하다.// SHIN (2006.2.4) e-15로 했을경우 오차범위가 너무 좁음
	if ( fabs(mathCross2D(vector1, vector2)) < dNormalZero*dTol && !( fabs(vector1[0] - vector2[0]) < dNormalZero*dTol && fabs(vector1[1] - vector2[1]) < dNormalZero*dTol ) )
	{
		return true;
	}

	return false;
}

bool CDgnMathFunc::mathIsInsidePoint2D(double p1[2], const int nData, double polyLine[][2], bool bIncludeOutLine = true)
{
	// p1 점에서 오른쪽(x+)으로 반직선을 그어서 폴리곤과 만나는 점이 짝수면 외보, 홀수면 내부
	double p2[2] = { p1[0]+100.0, p1[1] };

	for ( int i = 0; i < nData; i++ )
	{
		if ( p2[0] < polyLine[i][0] + 1.0e-8 )
		{
			p2[0] = polyLine[i][0] + 100.0; //p2의 X값은 polyline의 최대 X값보다 커야함
		}
	}

	// first check ; is On Line
	for ( int i = 0; i < nData; i++ )
	{
		int iIdx1 = i;
		int iIdx2 = i+1;
		if ( iIdx2 == nData ) iIdx2=0;

		if ( mathIsPointOfLine2D(polyLine[iIdx1], polyLine[iIdx2], p1) )
		{
			return bIncludeOutLine;
		}
	}

	// check polygon
	int nCount = 0;
	for ( int i = 0; i < nData; i++ )
	{
		int iIdx1 = i;
		int iIdx2 = i+1;
		if ( iIdx2 == nData ) iIdx2=0;

		int ic = mathIntersect_ccw2D(p1, p2, polyLine[iIdx1], polyLine[iIdx2]);
		if ( ic > 0 )
		{
			nCount++;
		}
		else if ( ic == 0 )
		{
			// check second point ; 폴리라인 양단 끝점에 걸리는 경우 반드시 2개의 폴리라인에 걸림. 따라서 1개만 count 하도록 2번째 point만 검토
			if ( IsNormalZero(p1[1] - __min(polyLine[iIdx1][1], polyLine[iIdx2][1])) && ( polyLine[iIdx1][1] != polyLine[iIdx2][1] ) )
			{
				nCount++;
			}
		}
	}

	return ( nCount % 2 ) == 1;
}

int CDgnMathFunc::mathIntersect_ccw2D(double p1Org[2], double p2Org[2], double p3Org[2], double p4Org[2])
{
	double p1[2], p2[2], p3[2], p4[2];
	for ( int i = 0; i < 2; i++ )
	{
		p1[i] = p1Org[i]; p2[i] = p2Org[i]; p3[i] = p3Org[i]; p4[i] = p4Org[i];
	}

	if ( p1[0] > p2[0] )
	{
		mathSwap2D(p1, p2);
	}
	if ( p3[0] > p4[0] )
	{
		mathSwap2D(p3, p4);
	}

	int r123 = math_ccw(p1[0], p1[1], p2[0], p2[1], p3[0], p3[1]);
	int r124 = math_ccw(p1[0], p1[1], p2[0], p2[1], p4[0], p4[1]);
	int r341 = math_ccw(p3[0], p3[1], p4[0], p4[1], p1[0], p1[1]);
	int r342 = math_ccw(p3[0], p3[1], p4[0], p4[1], p2[0], p2[1]);

	if ( r123*r124 < 0 && r341*r342 < 0 )
	{
		return 1; //교차하는 경우
	}

	if ( r123 == 0 && r124 == 0 )
	{
		if ( p1[0] != p2[0] )
		{
			if ( !( p3[0] > p2[0] || p1[0] > p4[0] ) )
			{
				return 0;
			}
			else
			{
				return -1;
			}
		}
		else
		{
			if ( !( __min(p3[1], p4[1]) > __max(p1[1], p2[1]) || __min(p1[1], p2[1]) > __max(p3[1], p4[1]) ) )
			{
				return 0;
			}
			else
			{
				return -1;
			}
		}
	}

	if ( r123 == 0 )
	{
		if ( p1[0] != p2[0] )
		{
			if ( p1[0] <= p3[0] && p3[0] <= p2[0] )
			{
				return 0;
			}
			else
			{
				return -1;
			}
		}
		else
		{
			if ( __min(p1[1], p2[1]) <= p3[1] && p3[1] <= __max(p1[1], p2[1]) )
			{
				return 0;
			}
			else
			{
				return -1;
			}
		}
	}
	if ( r124 == 0 )
	{
		if ( p1[0] != p2[0] )
		{
			if ( p1[0] <= p4[0] && p4[0] <= p2[0] )
			{
				return 0;
			}
			else
			{
				return -1;
			}
		}
		else
		{
			if ( __min(p1[1], p2[1]) <= p4[1] && p4[1] <= __max(p1[1], p2[1]) )
			{
				return 0;
			}
			else
			{
				return -1;
			}
		}
	}
	if ( r341 == 0 )
	{
		if ( p3[0] != p4[0] )
		{
			if ( p3[0] <= p1[0] && p1[0] <= p4[0] )
			{
				return 0;
			}
			else
			{
				return -1;
			}
		}
		else
		{
			if ( __min(p3[1], p4[1]) <= p1[1] && p1[1] <= __max(p3[1], p4[1]) )
			{
				return 0;
			}
			else
			{
				return -1;
			}
		}
	}
	if ( r342 == 0 )
	{
		if ( p3[0] != p4[0] )
		{
			if ( p3[0] <= p2[0] && p2[0] <= p4[0] )
			{
				return 0;
			}
			else
			{
				return -1;
			}
		}
		else
		{
			if ( __min(p3[1], p4[1]) <= p2[1] && p2[1] <= __max(p3[1], p4[1]) )
			{
				return 0;
			}
			else
			{
				return -1;
			}
		}
	}

	return -1;
}

int CDgnMathFunc::mathIntersect_ccw2D_Tol(double p1Org[2], double p2Org[2], double p3Org[2], double p4Org[2])
{
	double p1[2], p2[2], p3[2], p4[2];
	for ( int i = 0; i < 2; i++ )
	{
		p1[i] = p1Org[i]; p2[i] = p2Org[i]; p3[i] = p3Org[i]; p4[i] = p4Org[i];
	}

	if ( p1[0] > p2[0] )
	{
		mathSwap2D(p1, p2);
	}
	if ( p3[0] > p4[0] )
	{
		mathSwap2D(p3, p4);
	}

	int r123 = math_ccw_AutoTol(p1[0], p1[1], p2[0], p2[1], p3[0], p3[1]);
	int r124 = math_ccw_AutoTol(p1[0], p1[1], p2[0], p2[1], p4[0], p4[1]);
	int r341 = math_ccw_AutoTol(p3[0], p3[1], p4[0], p4[1], p1[0], p1[1]);
	int r342 = math_ccw_AutoTol(p3[0], p3[1], p4[0], p4[1], p2[0], p2[1]);

	if ( r123*r124 < 0 && r341*r342 < 0 )
		return 1; //교차하는 경우

	if ( r123 == 0 && r124 == 0 )
	{
		if ( fabs(p1[0] - p2[0]) > m_NormalZero )
		{
			if ( !( p3[0] > p2[0] || p1[0] > p4[0] ) ) return 0;
			else return -1;
		}
		else
		{
			if ( !( __min(p3[1], p4[1]) > __max(p1[1], p2[1]) || __min(p1[1], p2[1]) > __max(p3[1], p4[1]) ) ) return 0;
			else return -1;
		}
	}

	if ( r123 == 0 )
	{
		if ( fabs(p1[0] - p2[0]) > m_NormalZero ) 
		{
			if ( p1[0] < p3[0] + m_NormalZero && p3[0] < p2[0] + m_NormalZero ) return 0;
			else return -1;
		}
		else
		{
			if ( __min(p1[1], p2[1]) < p3[1] + m_NormalZero && p3[1] < __max(p1[1], p2[1]) + m_NormalZero ) return 0; 
			else return -1;
		}
	}
	if ( r124 == 0 )
	{
		if ( fabs(p1[0] - p2[0]) > m_NormalZero )
		{
			if ( p1[0] < p4[0] + m_NormalZero && p4[0] < p2[0] + m_NormalZero ) return 0;
			else return -1;
		}
		else
		{
			if ( __min(p1[1], p2[1]) < p4[1] + m_NormalZero && p4[1] < __max(p1[1], p2[1]) + m_NormalZero ) return 0;
			else return -1;
		}
	}
	if ( r341 == 0 )
	{
		if ( fabs(p3[0] - p4[0]) > m_NormalZero )
		{
			if ( p3[0] < p1[0] + m_NormalZero && p1[0] < p4[0] + m_NormalZero ) return 0;
			else return -1;
		}
		else
		{
			if ( __min(p3[1], p4[1]) < p1[1] + m_NormalZero && p1[1] < __max(p3[1], p4[1]) + m_NormalZero ) return 0;
			else return -1;
		}
	}
	if ( r342 == 0 )
	{
		if ( fabs(p3[0] - p4[0]) > m_NormalZero )
		{
			if ( p3[0] < p2[0] + m_NormalZero && p2[0] < p4[0] + m_NormalZero ) return 0;
			else return -1;
		}
		else
		{
			if ( __min(p3[1], p4[1]) < p2[1] + m_NormalZero && p2[1] < __max(p3[1], p4[1]) + m_NormalZero ) return 0;
			else return -1;
		}
	}

	return -1;
}

int CDgnMathFunc::math_ccw_AutoTol(double ax, double ay, double bx, double by, double cx, double cy)
{
	double dAbsMax = 0.;
	double bxcy = bx * cy; if ( fabs(dAbsMax) < fabs(bxcy) ) dAbsMax = bxcy;
	double aybx = ay * bx; if ( fabs(dAbsMax) < fabs(bxcy) ) dAbsMax = aybx;
	double axcy = ax * cy; if ( fabs(dAbsMax) < fabs(bxcy) ) dAbsMax = axcy;
	double bycx = by * cx; if ( fabs(dAbsMax) < fabs(bxcy) ) dAbsMax = bycx;
	double axby = ax * by; if ( fabs(dAbsMax) < fabs(bxcy) ) dAbsMax = axby;
	double aycx = ay * cx; if ( fabs(dAbsMax) < fabs(bxcy) ) dAbsMax = aycx;

	double l = bxcy - aybx - axcy - bycx + axby + aycx;

	// double의 곱연산을 고려한 유효자릿수로 Torerance 설정
	double dTol = fabs(dAbsMax * 1.e-9);
	if ( dTol == 0. ) dTol = m_NormalZero;

	if ( fabs(l) < dTol ) return 0;
	else if ( l > 0.0 ) return 1;
	else return -1;
}

// 직선과 직선의 교점 구하기
int CDgnMathFunc::mathLineLineCross2D(double line1[2][2], double line2[2][2], double cross[2])
{
	double vector1[2];
	vector1[0] = line1[0][0] - line1[1][0];
	vector1[1] = line1[0][1] - line1[1][1];
	mathNormalize2D(vector1, vector1);
	double vector2[2];
	vector2[0] = line2[0][0] - line2[1][0];
	vector2[1] = line2[0][1] - line2[1][1];
	mathNormalize2D(vector2, vector2);

	double product = mathCross2D(vector1, vector2);
	if ( fabs(product) < m_NormalZero * 100000.0 )   // 두 벡터가 평행인 경우
	{// SHIN (2006.2.4) e-15로 했을경우 오차범위가 너무 좁음
		cross[0] = 0.0;
		cross[1] = 0.0;
		return 0;
	}

	double vector3[2];
	vector3[0] = line1[0][0] - line2[0][0];
	vector3[1] = line1[0][1] - line2[0][1];
	double H = mathCross2D(vector2, vector3);
	cross[0] = line1[0][0] + vector1[0] * ( H / product );
	cross[1] = line1[0][1] + vector1[1] * ( H / product );
	return 1;
}

// 직선과 선분의 교점 구하기
int CDgnMathFunc::mathLineSegCross2D(double line[2][2], double sttP[2], double endP[2], double cross[2])
{
	double line2[2][2];
	line2[0][0] = sttP[0];
	line2[0][1] = sttP[1];
	line2[1][0] = endP[0];
	line2[1][1] = endP[1];

	if ( mathLineLineCross2D(line, line2, cross) == 0 )
	{
		return 0;
	}

	if ( mathIsPointOfLine2D(sttP, endP, cross) )
	{
		return 1;
	}
	return 0;
}

bool CDgnMathFunc::mathOffsetOfPolyline(double dOffset, const int nData, double polyLine[][2], int nRotType)
{
	// if ( nData < 3 ) return false;

	// double(*polyLine_Org)[2];
	// polyLine_Org = new double[nData][2];
	// int i, j;
	// for ( i = 0; i < nData; i++ )
	// {
	// 	for ( j = 0; j < 2; j++ )
	// 	{
	// 		polyLine_Org[i][j] = polyLine[i][j];
	// 	}
	// }
	// if ( nRotType <= 0 )
	// {// Polyline의 회전방향을 자동으로 검토
	// 	nRotType = GetRotationType(nData, polyLine_Org);
	// }

	// double dZero = m_NormalZero * 1000.0; //동일점 여부 판별시 오차법위를 넓혀줌
	// double vector[2];
	// double cross[2];
	// double vectorOuter[2];//바깥방향으로의 단위벡터
	// double line1[2][2];
	// double line2[2][2];
	// double dL_line1, dL_line2;
	// int nStrID, nEndID;

	// //[nData-1] ~ [0] 선분
	// nStrID = nData - 1;
	// nEndID = 0;
	// for ( i = nData - 1; i > 0; i-- )
	// {
	// 	dL_line1 = mathLength(polyLine_Org[nEndID][0] - polyLine_Org[i][0], polyLine_Org[nEndID][1] - polyLine_Org[i][1]);
	// 	if ( dL_line1 >= m_NormalZero )
	// 	{
	// 		nStrID = i;
	// 		break;
	// 	}
	// }

	// vector[0] = polyLine_Org[nEndID][0] - polyLine_Org[nStrID][0];
	// vector[1] = polyLine_Org[nEndID][1] - polyLine_Org[nStrID][1];
	// mathNormalize(vector[0], vector[1], vector[0], vector[1]);
	// if ( nRotType == 1 )
	// {
	// 	vectorOuter[0] = vector[1];
	// 	vectorOuter[1] = -1.0*vector[0];
	// }
	// else
	// {
	// 	vectorOuter[0] = -1.0*vector[1];
	// 	vectorOuter[1] = vector[0];
	// }
	// line1[0][0] = polyLine_Org[nStrID][0] + dOffset * vectorOuter[0];
	// line1[0][1] = polyLine_Org[nStrID][1] + dOffset * vectorOuter[1];
	// line1[1][0] = polyLine_Org[nEndID][0] + dOffset * vectorOuter[0];
	// line1[1][1] = polyLine_Org[nEndID][1] + dOffset * vectorOuter[1];
	// dL_line1 = mathLength(line1[0][0], line1[0][1], line1[1][0], line1[1][1]);
	// CArray<int, int> arSameID;
	// int nSameSize = 0;
	// for ( i = 0; i < nData; i++ )
	// {
	// 	//[i] ~ [i+1] 선분
	// 	if ( i < nData - 1 )
	// 	{
	// 		nStrID = i;
	// 		nEndID = i + 1;
	// 	}
	// 	else
	// 	{
	// 		nStrID = i;
	// 		nEndID = 0;
	// 	}

	// 	vector[0] = polyLine_Org[nEndID][0] - polyLine_Org[nStrID][0];
	// 	vector[1] = polyLine_Org[nEndID][1] - polyLine_Org[nStrID][1];
	// 	mathNormalize(vector[0], vector[1], vector[0], vector[1]);
	// 	if ( nRotType == 1 )
	// 	{
	// 		vectorOuter[0] = vector[1];
	// 		vectorOuter[1] = -1.0*vector[0];
	// 	}
	// 	else
	// 	{
	// 		vectorOuter[0] = -1.0*vector[1];
	// 		vectorOuter[1] = vector[0];
	// 	}
	// 	line2[0][0] = polyLine_Org[nStrID][0] + dOffset * vectorOuter[0];
	// 	line2[0][1] = polyLine_Org[nStrID][1] + dOffset * vectorOuter[1];
	// 	line2[1][0] = polyLine_Org[nEndID][0] + dOffset * vectorOuter[0];
	// 	line2[1][1] = polyLine_Org[nEndID][1] + dOffset * vectorOuter[1];
	// 	dL_line2 = mathLength(line2[0][0], line2[0][1], line2[1][0], line2[1][1]);
	// 	if ( dL_line2 < m_NormalZero )
	// 	{
	// 		arSameID.Add(i);
	// 		if ( i == nData - 1 )
	// 		{
	// 			nSameSize = arSameID.GetSize();
	// 			for ( j = 0; j < nSameSize; j++ )
	// 			{
	// 				polyLine[arSameID[j]][0] = polyLine[0][0];
	// 				polyLine[arSameID[j]][1] = polyLine[0][1];
	// 			}
	// 			arSameID.RemoveAll();
	// 		}
	// 	}
	// 	else
	// 	{
	// 		if ( mathLength(line1[1][0], line1[1][1], line2[0][0], line2[0][1]) <= dZero * ( dL_line1 + dL_line2 ) )
	// 		{// 직선상에 있을 경우에는 line1의 끝점과 line2의 시작점이 동일
	// 			polyLine[i][0] = line2[0][0];
	// 			polyLine[i][1] = line2[0][1];
	// 		}
	// 		else
	// 		{
	// 			if ( mathLineLineCross2D(line1, line2, cross) == 0 )
	// 			{//직선상이 아닌 평행일 경우 계산불가
	// 				delete[] polyLine_Org;
	// 				return false;
	// 			}
	// 			else
	// 			{
	// 				polyLine[i][0] = cross[0];
	// 				polyLine[i][1] = cross[1];
	// 			}
	// 		}

	// 		if ( arSameID.GetSize() > 0 )
	// 		{
	// 			nSameSize = arSameID.GetSize();
	// 			for ( j = 0; j < nSameSize; j++ )
	// 			{
	// 				polyLine[arSameID[j]][0] = polyLine[i][0];
	// 				polyLine[arSameID[j]][1] = polyLine[i][1];
	// 			}
	// 			arSameID.RemoveAll();
	// 		}

	// 		// line2정보를 다음 좌표점계산을 위해 line1정보에 대입
	// 		line1[0][0] = line2[0][0];
	// 		line1[0][1] = line2[0][1];
	// 		line1[1][0] = line2[1][0];
	// 		line1[1][1] = line2[1][1];
	// 		dL_line1 = dL_line2;
	// 	}
	// }
	// delete[] polyLine_Org;
	return true;
}

int CDgnMathFunc::mathCircleForm3Point(double p1[2], double p2[2], double p3[2], double CenterP[2], double& dRadius)
{
	CenterP[0] = CenterP[1] = dRadius = 0.0;
	double dZero = m_NormalZero * 100000.0;
	double dDis12 = mathLength(p1[0], p1[1], p2[0], p2[1]);
	double dDis13 = mathLength(p1[0], p1[1], p3[0], p3[1]);
	double dDis23 = mathLength(p2[0], p2[1], p3[0], p3[1]);

	if ( dDis12 < dZero && dDis13 < dZero ) return -1;
	if ( dDis12 < dZero || dDis13 < dZero || dDis23 < dZero ) return 0;

	double vector1[2];
	vector1[0] = p1[0] - p2[0];
	vector1[1] = p1[1] - p2[1];
	mathNormalize2D(vector1, vector1);
	double vector2[2];
	vector2[0] = p3[0] - p2[0];
	vector2[1] = p3[1] - p2[1];
	mathNormalize2D(vector2, vector2);

	double product = mathCross2D(vector1, vector2);
	if ( fabs(product) < dZero ) return 0;  // 두 벡터가 평행인 경우

	// 1,2점의 중심점에서 수직인 선과 2,3점의 중심점에서 수직인 선과의 교차점이 원의 중심점임
	double line1[2][2], line2[2][2];
	line1[0][0] = ( p1[0] + p2[0] ) / 2.0;
	line1[0][1] = ( p1[1] + p2[1] ) / 2.0;
	line1[1][0] = line1[0][0] + vector1[1];
	line1[1][1] = line1[0][1] - vector1[0];
	line2[0][0] = ( p3[0] + p2[0] ) / 2.0;
	line2[0][1] = ( p3[1] + p2[1] ) / 2.0;
	line2[1][0] = line2[0][0] + vector2[1];
	line2[1][1] = line2[0][1] - vector2[0];

	if ( mathLineLineCross2D(line1, line2, CenterP) != 1 ) return 0;  // 두 벡터가 평행인 경우

	dRadius = mathLength(p1[0], p1[1], CenterP[0], CenterP[1]);
	return 1;
}

int CDgnMathFunc::mathArcForm3Point(double p1[2], double p2[2], double p3[2], double CenterP[2], double& dRadius, double& dStartAngle, double& dSweepAngle)
{
	dStartAngle = dSweepAngle = 0.0;
	int nChk = mathCircleForm3Point(p1, p2, p3, CenterP, dRadius);
	if ( nChk != 1 ) return nChk;

	double vector0[2];
	vector0[0] = 1.0;
	vector0[1] = 0.0;
	double vector1[2];
	vector1[0] = p1[0] - CenterP[0];
	vector1[1] = p1[1] - CenterP[1];
	mathNormalize2D(vector1, vector1);
	double vector2[2];
	vector2[0] = p2[0] - CenterP[0];
	vector2[1] = p2[1] - CenterP[1];
	mathNormalize2D(vector2, vector2);
	double vector3[2];
	vector3[0] = p3[0] - CenterP[0];
	vector3[1] = p3[1] - CenterP[1];
	mathNormalize2D(vector3, vector3);

	double dAngle1 = mathCrossAngle2DSign(vector0[0], vector0[1], vector1[0], vector1[1]);
	double dAngle2 = mathCrossAngle2DSign(vector1[0], vector1[1], vector2[0], vector2[1]);
	double dAngle3 = mathCrossAngle2DSign(vector1[0], vector1[1], vector3[0], vector3[1]);

	dStartAngle = dAngle1;
	if ( dAngle2*dAngle3 >= 0.0 )
	{
		if ( fabs(dAngle3) >= fabs(dAngle2) )
		{
			dSweepAngle = dAngle3;
		}
		else
		{
			if ( dAngle3 >= 0.0 )
			{
				dSweepAngle = dAngle3 - 360.0;
			}
			else
			{
				dSweepAngle = dAngle3 + 360.0;
			}
		}
	}
	else
	{
		if ( dAngle3 >= 0.0 )
		{
			dSweepAngle = dAngle3 - 360.0;
		}
		else
		{
			dSweepAngle = dAngle3 + 360.0;
		}
	}
	return 1;
}

int CDgnMathFunc::mathArcForm3Point(double p1[2], double p2[2], double p3[2], double CenterP[2], double& dRadius, double& dStartAngle, double& dSweepAngle, double& dBulge)
{
	dBulge = 0.0;
	int nChk = mathArcForm3Point(p1, p2, p3, CenterP, dRadius, dStartAngle, dSweepAngle);
	dBulge = dRadius * ( 1.0 - cos(fabs(dSweepAngle)*m_trrad) );
	if ( dSweepAngle < 0.0 )
	{
		dBulge *= ( -1.0 );
	}
	return nChk;
}

bool CDgnMathFunc::GetCenter_Bulge(double line[2][2], double bulge, double centerP[2])
{
	centerP[0] = ( line[0][0] + line[1][0] ) / 2.0;
	centerP[1] = ( line[0][1] + line[1][1] ) / 2.0;
	if ( bulge == 0.0 ) return false;

	double h = mathLength(line[0][0], line[0][1], line[1][0], line[1][1]) / 2.0;
	double len = ( h*h - bulge * bulge ) / ( 2.0*bulge );

	double vector[2];
	mathNormalize(( line[0][0] - line[1][0] ), ( line[0][1] - line[1][1] ), vector[0], vector[1]);
	// 단위백터를 90회전된 각도로 계산한다.
	centerP[0] += len * vector[1];
	centerP[1] -= len * vector[0];
	return true;
}

double CDgnMathFunc::GetBulge(double line[2][2], double radius, bool inside)
{
	int sign = radius >= 0.0 ? 1 : -1;
	double absR = fabs(radius);

	double h = mathLength(line[0][0], line[0][1], line[1][0], line[1][1]) / 2.0;
	if ( h > absR + m_NormalZero )
	{
		return 0;
	}
	else if ( IsNormalZero(h - absR) )
	{
		return radius;
	}

	double bulge = ( absR - sqrt(absR*absR - h * h) ) * sign;
	if ( !inside )
	{
		bulge = radius * 2.0 - bulge;
	}
	return bulge;
}

int CDgnMathFunc::GetRotationType(const int nData, double polyLine[][2])
{
	int nRotType = 0;
	//진행방향에서 좌.우중 안쪽방향을 알기위해
	// 첫번째 선분의 중심점에서 진행방향에서 좌측편으로 약간 떨어진 점이 polyline에 포함되는지 여부로 검토
	double vector[2];
	vector[0] = polyLine[1][0] - polyLine[0][0];
	vector[1] = polyLine[1][1] - polyLine[0][1];
	mathNormalize(vector[0], vector[1], vector[0], vector[1]);
	double vectorLeft[2];
	vectorLeft[0] = -1.0*vector[1];
	vectorLeft[1] = vector[0];
	double dDeltaLength = 1.e-6 * mathLength(polyLine[1][0] - polyLine[0][0], polyLine[1][1] - polyLine[0][1]);
	double dPoint[2];
	dPoint[0] = ( polyLine[1][0] + polyLine[0][0] ) / 2.0 + dDeltaLength * vectorLeft[0];
	dPoint[1] = ( polyLine[1][1] + polyLine[0][1] ) / 2.0 + dDeltaLength * vectorLeft[1];
	bool bInside = mathIsInsidePoint2D(dPoint, nData, polyLine, true);
	if ( bInside )
	{
		nRotType = 1;
	}
	else
	{
		nRotType = 2;
	}

	if ( nRotType == 0 );

	return nRotType;
}

bool CDgnMathFunc::Get_VolmCenter(double pointCenV[3], double& dVolm, const int nData, double PontUnitList[][3], int nAxis)
{
	// ※주의사항 : 기둥형상에 쐐기 또는 갈라진 틈과 같은 모양이 있을 경우 (폭-길이비 1/1000000000000000 이하)
	//              단면2차모멘트 값이 잘못된 값을 나타낼 수 있음

	pointCenV[0] = pointCenV[1] = pointCenV[2] = dVolm = 0.0;

	if ( nData < 3 ) return false;

	int ix = 0, iy = 0, iz = 0;
	switch ( nAxis )
	{
	case 1: { ix = 0; iy = 1; iz = 2; } break; // xy평면
	case 2: { ix = 2; iy = 0; iz = 1; } break; // yz평면
	case 3: { ix = 0; iy = 2; iz = 1; } break; // xz평면
	}

	// x가 0인 y-z평면에서 부터의 arPontUnit좌표까지의 기둥에 대한 부피 및 부피의 도심점을 산출하며
	// 이때 arPontUnit좌표값들은 하나의 평면상에 존재 하여야 함 (즉 해당 기둥을 수평 또는 경사지게 자른 모양이어야 함)
	// [계산방법]
	// 시작점을 포함하는 연결된 3각기둥들로 나누어서 각각의 체적및 도심들을 구하여 전체의 값을 구함
	double dpCenV[3];

	int iNodeNum = nData;
	int i = 0, j = 0, k = 0;

	// 계산할 좌표순서는 반시계 방향으로 정의 되어 있어야 함=============================================
	// Polyline의 회전방향을 검토
	// 진행방향에서 좌.우중 안쪽방향을 알기위해
	// 첫번째 선분의 중심점에서 진행방향에서 좌측편으로 약간 떨어진 점이 polyline에 포함되는지 여부로 검토
	double(*polyLine)[2];
	polyLine = new double[iNodeNum][2];
	int iNode_MaxLeng = 0;
	double dMaxLeng = 0.0;
	double dLength;
	for ( i = 0; i < iNodeNum; i++ )
	{
		polyLine[i][0] = PontUnitList[i][ix];
		polyLine[i][1] = PontUnitList[i][iy];

		if ( i == iNodeNum - 1 )
		{
			dLength = CDgnMathFunc::mathLength(PontUnitList[i][ix] - PontUnitList[0][ix], PontUnitList[i][iy] - PontUnitList[0][iy]);
		}
		else
		{
			dLength = CDgnMathFunc::mathLength(PontUnitList[i][ix] - PontUnitList[i + 1][ix], PontUnitList[i][iy] - PontUnitList[i + 1][iy]);
		}

		if ( dLength > dMaxLeng )
		{
			iNode_MaxLeng = i;
		}
	}

	double vector[2], pointA[2], pointB[2];
	pointA[0] = polyLine[iNode_MaxLeng][0];
	pointA[1] = polyLine[iNode_MaxLeng][1];
	if ( iNode_MaxLeng == iNodeNum - 1 )
	{
		pointB[0] = polyLine[0][0];
		pointB[1] = polyLine[0][1];
	}
	else
	{
		pointB[0] = polyLine[iNode_MaxLeng + 1][0];
		pointB[1] = polyLine[iNode_MaxLeng + 1][1];
	}
	vector[0] = pointB[0] - pointA[0];
	vector[1] = pointB[1] - pointA[1];
	CDgnMathFunc::mathNormalize(vector[0], vector[1], vector[0], vector[1]);
	double vectorLeft[2];
	vectorLeft[0] = -1.0*vector[1];
	vectorLeft[1] = vector[0];
	double dDeltaLength = 1.e-6 * CDgnMathFunc::mathLength(pointB[0] - pointA[0], pointB[1] - pointA[1]);
	double dPoint[2];
	dPoint[0] = ( pointB[0] + pointA[0] ) / 2.0 + dDeltaLength * vectorLeft[0];
	dPoint[1] = ( pointB[1] + pointA[1] ) / 2.0 + dDeltaLength * vectorLeft[1];
	bool bInside = CDgnMathFunc::mathIsInsidePoint2D(dPoint, iNodeNum, polyLine, true);
	if ( polyLine ) delete[] polyLine;

	//반시계 방향으로 좌표 Sort
	double(*PUnitList_Sort)[3];
	PUnitList_Sort = new double[iNodeNum][3];
	if ( bInside )
	{//반시계 방향
		for ( i = 0; i < iNodeNum; i++ )
		{
			for ( j = 0; j < 3; j++ )
			{
				PUnitList_Sort[i][j] = PontUnitList[i][j];
			}
		}
	}
	else
	{//시계 방향
		for ( i = 0; i < iNodeNum; i++ )
		{
			for ( j = 0; j < 3; j++ )
			{
				PUnitList_Sort[i][j] = PontUnitList[iNodeNum - 1 - i][j];
			}
		}
	}
	//=====================================================================================================

	// dpCenV (Center of Volume).
	dVolm = 0.0;

	int iIterNum = iNodeNum - 2;

	double  *adTriVol;
	double(*adTriCen)[3];
	adTriVol = new double[iIterNum];
	adTriCen = new double[iIterNum][3];

	double dTriCen[3][3] =
	{
		{ 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0 }
	};

	// Get Volume.
	//CArray<_DB_PONT_UNIT,_DB_PONT_UNIT> arTriPontUnit;
	for ( i = 0; i < iIterNum; i++ )
	{
		double dWeight = 1.0 / 3.0;
		double dpFunc[3] = { 1.0, 1.0, 1.0 };
		double dpTransFunc[3][3] = { { 0.50, 0.25, 0.25 }, { 0.25, 0.50, 0.25 }, { 0.25, 0.25, 0.50 } };
		double dTriArea = 0.0, dTriVolm = 0.0;
		double dpTriCenV[3] = { 0.0, 0.0, 0.0 };

		// Get Volume by Triangual.
		double dpNode[3][3] =
		{
			{ PUnitList_Sort[0][ix], PUnitList_Sort[0][iy], PUnitList_Sort[0][iz] },
			{ PUnitList_Sort[i + 1][ix], PUnitList_Sort[i + 1][iy], PUnitList_Sort[i + 1][iz] },
			{ PUnitList_Sort[i + 2][ix], PUnitList_Sort[i + 2][iy], PUnitList_Sort[i + 2][iz] }
		};
		dTriArea = ( ( dpNode[1][0] - dpNode[0][0] ) * ( dpNode[2][1] - dpNode[0][1] ) - ( dpNode[1][1] - dpNode[0][1] ) * ( dpNode[2][0] - dpNode[0][0] ) ) / 2.0;// 3점이 반시계방향이면 '+', 시계방향이면  '-'값을 가짐

		if ( IsNormalZero(dTriArea) )
		{
			adTriVol[i] = 0.0;
			adTriCen[i][0] = adTriCen[i][1] = adTriCen[i][2] = 0.0;
			continue;
		}
		// Volume.
		dWeight *= dTriArea;
		for ( k = 0; k < 3; k++ )
		{
			dpFunc[k] *= dWeight;
		}
		for ( k = 0; k < 3; k++ )
		{
			dTriVolm += dpFunc[k] * dpNode[k][2];
		}
		// Centroid of Mass.
		for ( k = 0; k < 3; k++ )
		{
			for ( j = 0; j < 3; j++ )
			{
				dpTransFunc[k][j] *= dWeight;
			}
		}

		double dTemp[3] = { 0.0, 0.0, 0.0 };
		for ( k = 0; k < 3; k++ )
		{
			for ( j = 0; j < 3; j++ )
			{
				dTemp[k] += dpTransFunc[k][j] * dpNode[j][2];
			}
		}
		for ( k = 0; k < 3; k++ )
		{
			dpTriCenV[k] = 0.0;
			for ( j = 0; j < 3; j++ )
			{
				dpTriCenV[k] += dpNode[j][k] * dTemp[j];
			}
			dpTriCenV[k] /= dTriVolm;
		}
		dpTriCenV[2] /= 2.0;

		adTriVol[i] = dTriVolm;
		adTriCen[i][0] = dpTriCenV[0];
		adTriCen[i][1] = dpTriCenV[1];
		adTriCen[i][2] = dpTriCenV[2];

		dVolm += dTriVolm;
	}
	// Get Centroid.
	for ( i = 0; i < 3; i++ )
	{
		dpCenV[i] = 0.0;
		for ( j = 0; j < iIterNum; j++ )
		{
			dpCenV[i] += adTriCen[j][i] * adTriVol[j];
		}
		dpCenV[i] /= dVolm;
	}

	if ( adTriVol ) delete[]adTriVol;
	if ( adTriCen ) delete[]adTriCen;
	if ( PUnitList_Sort ) delete[]PUnitList_Sort;

	pointCenV[ix] = dpCenV[0];
	pointCenV[iy] = dpCenV[1];
	pointCenV[iz] = dpCenV[2];
	return true;
}

bool CDgnMathFunc::IsOverLap_Polyline(const int nData1, double polyLine1[][2], const int nData2, double polyLine2[][2])
{
	if ( nData1 == 0 || nData2 == 0 ) return false;

	int i, j;
	// 교차 검토
	double p1Org[2], p2Org[2], p3Org[2], p4Org[2];
	for ( i = 0; i < nData1; i++ )
	{
		p1Org[0] = polyLine1[i][0];
		p1Org[1] = polyLine1[i][1];
		p2Org[0] = polyLine1[( i == nData1 - 1 ? 0 : i )][0];
		p2Org[1] = polyLine1[( i == nData1 - 1 ? 0 : i )][1];

		for ( j = 0; j < nData2; j++ )
		{
			p3Org[0] = polyLine2[j][0];
			p3Org[1] = polyLine2[j][1];
			p4Org[0] = polyLine2[( j == nData2 - 1 ? 0 : j )][0];
			p4Org[1] = polyLine2[( j == nData2 - 1 ? 0 : j )][1];

			if ( mathIntersect_ccw2D(p1Org, p2Org, p3Org, p4Org) >= 0 ) return true; //교차하거나 포함되는 경우
		}
	}

	// 내부점 검토
	double dxmax1, dxmin1, dymax1, dymin1;
	double dxmax2, dxmin2, dymax2, dymin2;
	for ( i = 0; i < nData1; i++ )
	{
		if ( i == 0 )
		{
			dxmax1 = dxmin1 = polyLine1[i][0];
			dymax1 = dymin1 = polyLine1[i][1];
		}
		else
		{
			dxmax1 = __max(dxmax1, polyLine1[i][0]);
			dxmin1 = __min(dxmin1, polyLine1[i][0]);
			dymax1 = __max(dymax1, polyLine1[i][1]);
			dymin1 = __min(dymin1, polyLine1[i][1]);
		}
	}
	for ( i = 0; i < nData2; i++ )
	{
		if ( i == 0 )
		{
			dxmax2 = dxmin2 = polyLine2[i][0];
			dymax2 = dymin2 = polyLine2[i][1];
		}
		else
		{
			dxmax2 = __max(dxmax2, polyLine2[i][0]);
			dxmin2 = __min(dxmin2, polyLine2[i][0]);
			dymax2 = __max(dymax2, polyLine2[i][1]);
			dymin2 = __min(dymin2, polyLine2[i][1]);
		}
	}

	double p1[2];
	if ( __max(dxmax1 - dxmin1, dymax1 - dymin1) > __max(dxmax2 - dxmin2, dymax2 - dymin2) || ( nData1 > 1 && nData2 == 1 ) )
	{
		for ( i = 0; i < nData2; i++ )
		{
			p1[0] = polyLine2[i][0];
			p1[1] = polyLine2[i][1];
			if ( mathIsInsidePoint2D(p1, nData1, polyLine1, true) ) return true;
		}
	}
	else if ( nData2 > 1 )
	{
		for ( i = 0; i < nData1; i++ )
		{
			p1[0] = polyLine1[i][0];
			p1[1] = polyLine1[i][1];
			if ( mathIsInsidePoint2D(p1, nData2, polyLine2, true) ) return true;
		}
	}
	return false;
}

double CDgnMathFunc::mathSin(double x)
{
	double r = sin(x);
	if ( IsNormalZero(r) ) r =  0.0;
	if ( IsNormalZero(1.0 - r) ) r =  1.0;
	if ( IsNormalZero(1.0 + r) ) r = -1.0;
	return r;
}

double CDgnMathFunc::mathCos(double x)
{
	double r = cos(x);
	if ( IsNormalZero(r) ) r =  0.0;
	if ( IsNormalZero(1.0 - r) ) r =  1.0;
	if ( IsNormalZero(1.0 + r) ) r = -1.0;
	return r;
}

double CDgnMathFunc::mathTan(double x)
{
	double dsin = mathSin(x);
	double dcos = mathCos(x);
	if ( IsNormalZero(dcos) ) return DBL_MAX;
	return dsin / dcos;
}

double CDgnMathFunc::mathCot(double x)
{
	double dtan = mathTan(x);
	if ( IsNormalZero(dtan) ) return DBL_MAX;
	return 1.0 / dtan;
}

double CDgnMathFunc::mathAsin(double x)
{
	if ( x < -1 ) x = -1.0;
	if ( x >  1 ) x =  1.0;
	return asin(x);
}

double CDgnMathFunc::mathAcos(double x)
{
	if ( x < -1 ) x = -1.0;
	if ( x >  1 ) x =  1.0;
	return acos(x);
}

double CDgnMathFunc::mathSqrt(double x)
{
	if ( x < 0 ) x = 0.0;
	return sqrt(x);
}

int CDgnMathFunc::mathAbs(int x)
{
	return abs(x);
}

double CDgnMathFunc::mathAbs(double x)
{
	return fabs(x);
}

int CDgnMathFunc::RealToInt(const double dVal, double dTol)
{
	if ( dTol <= 0 );

	int nInt = ( int ) dVal;
	if ( fabs(dVal - nInt - 1) < dTol )
	{
		nInt = ( int ) ceil(dVal);
	}
	return nInt;
}

double CDgnMathFunc::mathMomentOfInertiaRectangle(double dB, double dH, double de)
{
	return dB * pow(dH, 3.0) / 12.0 + dB * dH*pow(de, 2.0);
}

double CDgnMathFunc::mathMomentOfInertiaTriangle(double dB, double dH, double de)
{
	return dB * pow(dH, 3.0) / 36.0 + 0.5*dB*dH*pow(de, 2.0);
}

double CDgnMathFunc::mathMomentOfInertiaCircle(double dr)
{
	return CDgnMathFunc::Get_pi()*pow(dr, 4.0) / 4.0;
}

double CDgnMathFunc::mathRoundOff(double dVal, int N)
{
	dVal *= pow(10, N);
	dVal += 0.5;
	return floor(dVal) / pow(10, N);
}

double CDgnMathFunc::mathRoundUp(double dVal, int N)
{
	dVal *= pow(10, N);
	dVal += 1.0;
	return floor(dVal) / pow(10, N);
}

double CDgnMathFunc::mathRoundDown(double dVal, int N)
{
	if ( N < 0 )
	{
		N *= -1;
		dVal /= pow(10, N);
		return floor(dVal) * pow(10, N);
	}
	else
	{
		dVal *= pow(10, N);
		return floor(dVal) / pow(10, N);
	}
	return 0.0;
}

bool CDgnMathFunc::mathFrameLocalVector(double Coor_i[3], double Coor_j[3], double Angle, double FrameLocalVector[3][3])
{
	double dx = Coor_j[0] - Coor_i[0];
	double dy = Coor_j[1] - Coor_i[1];
	double dz = Coor_j[2] - Coor_i[2];
	double Length = CDgnMathFunc::mathLength(dx, dy, dz);
	if ( Length < CDgnMathFunc::Get_NormalZero() ) return false;

	double Distance = CDgnMathFunc::mathLength(dx, dy);
	bool bVertical = ( Distance/Length < m_VertTol ) ? true : false;

	double xyz[3][3] = { 0.0, };
	xyz[2][0] = Coor_i[0];
	xyz[2][1] = Coor_i[1];
	xyz[2][2] = Coor_i[2];
	if ( bVertical )
	{
		xyz[2][0] += Length;
	}
	else
	{
		xyz[2][2] += Length;
	}
	if ( fabs(Angle) > m_AngleZero )
	{
		double ux = dx;
		double uy = dy;
		double uz = dz;
		CDgnMathFunc::mathNormalize(ux, uy, uz, ux, uy, uz);
		double x = xyz[2][0] - Coor_i[0];
		double y = xyz[2][1] - Coor_i[1];
		double z = xyz[2][2] - Coor_i[2];
		CDgnMathFunc::mathRotate(Angle, ux, uy, uz, x, y, z);
		xyz[2][0] = x + Coor_i[0];
		xyz[2][1] = y + Coor_i[1];
		xyz[2][2] = z + Coor_i[2];
	}

	xyz[0][0] = Coor_j[0] - Coor_i[0];
	xyz[0][1] = Coor_j[1] - Coor_i[1];
	xyz[0][2] = Coor_j[2] - Coor_i[2];
	xyz[1][0] = xyz[2][0] - Coor_i[0];
	xyz[1][1] = xyz[2][1] - Coor_i[1];
	xyz[1][2] = xyz[2][2] - Coor_i[2];

	if ( !CDgnMathFunc::mathNormalize(xyz[0], FrameLocalVector[0]) )return false;

	if ( !CDgnMathFunc::mathNormalize(xyz[1], FrameLocalVector[2]) )return false;

	CDgnMathFunc::mathNormal(FrameLocalVector[2], FrameLocalVector[0], FrameLocalVector[1]);
	if ( !CDgnMathFunc::mathNormalize(FrameLocalVector[1], FrameLocalVector[1]) )return false;

	CDgnMathFunc::mathNormal(FrameLocalVector[0], FrameLocalVector[1], FrameLocalVector[2]);
	return true;
}

bool CDgnMathFunc::mathPlaneLocalVector(int NumNode, double Coor[][3], double dAngle, double PlaneLocalVector[3][3])
{
	if ( NumNode == 3 )
	{
		PlaneLocalVector[0][0] = Coor[1][0] - Coor[0][0];
		PlaneLocalVector[0][1] = Coor[1][1] - Coor[0][1];
		PlaneLocalVector[0][2] = Coor[1][2] - Coor[0][2];
		PlaneLocalVector[1][0] = Coor[2][0] - Coor[1][0];
		PlaneLocalVector[1][1] = Coor[2][1] - Coor[1][1];
		PlaneLocalVector[1][2] = Coor[2][2] - Coor[1][2];
	}
	else if ( NumNode >= 4 )
	{
		PlaneLocalVector[0][0] = Coor[1][0] + Coor[2][0] - Coor[3][0] - Coor[0][0];
		PlaneLocalVector[0][1] = Coor[1][1] + Coor[2][1] - Coor[3][1] - Coor[0][1];
		PlaneLocalVector[0][2] = Coor[1][2] + Coor[2][2] - Coor[3][2] - Coor[0][2];
		PlaneLocalVector[1][0] = Coor[2][0] + Coor[3][0] - Coor[0][0] - Coor[1][0];
		PlaneLocalVector[1][1] = Coor[2][1] + Coor[3][1] - Coor[0][1] - Coor[1][1];
		PlaneLocalVector[1][2] = Coor[2][2] + Coor[3][2] - Coor[0][2] - Coor[1][2];
	}
	else
	{
		return false;
	}

	if ( !CDgnMathFunc::mathNormalize(PlaneLocalVector[0], PlaneLocalVector[0]) )return false;
	if ( !CDgnMathFunc::mathNormalize(PlaneLocalVector[1], PlaneLocalVector[1]) )return false;
	CDgnMathFunc::mathNormal(PlaneLocalVector[0], PlaneLocalVector[1], PlaneLocalVector[2]);
	if ( !CDgnMathFunc::mathNormalize(PlaneLocalVector[2], PlaneLocalVector[2]) )return false;
	CDgnMathFunc::mathNormal(PlaneLocalVector[2], PlaneLocalVector[0], PlaneLocalVector[1]);

	if ( fabs(dAngle) > CDgnMathFunc::m_AngleZero )
	{
		CDgnMathFunc::mathRotate(dAngle, 0., 0., 0., PlaneLocalVector[2][0], PlaneLocalVector[2][1], PlaneLocalVector[2][2], PlaneLocalVector[0][0], PlaneLocalVector[0][1], PlaneLocalVector[0][2]);
		CDgnMathFunc::mathRotate(dAngle, 0., 0., 0., PlaneLocalVector[2][0], PlaneLocalVector[2][1], PlaneLocalVector[2][2], PlaneLocalVector[1][0], PlaneLocalVector[1][1], PlaneLocalVector[1][2]);
	}

	return true;
}

bool CDgnMathFunc::IntersectTriangle(double RayOrigin[3], double RayDirection[3], double V0[3], double V1[3], double V2[3], double cross[3])
{
	dgn::math::CVector3 vRayOrigin(RayOrigin);
	dgn::math::CVector3 vRayDirection(RayDirection);
	dgn::math::CVector3 vV0(V0);
	dgn::math::CVector3 vV1(V1);
	dgn::math::CVector3 vV2(V2);
	dgn::math::CVector3 vCross(cross);

	// 1. make line V0-V1, V0-V2
	double edge1[3] = { 0.0, }, edge2[3] = { 0.0, };
	for ( int i = 0; i < 3; ++i )
	{
		edge1[i] = V1[i] - V0[i];
		edge2[i] = V2[i] - V0[i];
	}

	double pvec[3] = { 0.0, };
	CDgnMathFunc::mathCross(RayOrigin, edge2, pvec);
	double dot = CDgnMathFunc::mathDot(edge1, pvec, false);

	double tvec[3];
	if ( dot > 0.0 )
	{
		for ( int i=0; i<3; ++i )
			tvec[i] = RayOrigin[i] - V0[i];
	}
	else
	{
		for ( int i = 0; i < 3; ++i )
			tvec[i] = V0[i] - RayOrigin[i];
		dot = -dot;
	}

	if ( dot < 0.0001 )
		return false;

	double u = CDgnMathFunc::mathDot(tvec, pvec, false);
	// check finite plane.
	if ( u < 0.0 || u > dot )
		return false;

	double qvec[3] = { 0.0, };
	CDgnMathFunc::mathCross(tvec, edge1, qvec);

	double v = CDgnMathFunc::mathDot(RayDirection, qvec, false);
	// check finite plane.
	if ( v < 0.0 || ( u + v ) > dot ) return false;

	double t = CDgnMathFunc::mathDot(edge2, qvec, false);

	double flnvDet = 1.0 / dot;
	t *= flnvDet;
	u *= flnvDet;
	v *= flnvDet;

	// calc cross point
	for ( int i = 0; i < 3; ++i )
	{
		cross[i] = RayOrigin[i] + RayDirection[i] * t;
	}

	return true;
}

bool CDgnMathFunc::IntersectTriangle(const CVector3& RayOrigin, const CVector3& RayDirection, const CVector3& V0, const CVector3& V1, const CVector3& V2, CVector3& cross)
{
	CVector3 edge1(V1 - V0), edge2(V2 - V0);
	CVector3 pvec(RayOrigin);
	pvec.cross(pvec);
	double dot = edge1.dot(pvec);

	CVector3 tvec;
	if ( dot > 0.0 )
	{
		tvec = RayOrigin - V0;
	}
	else
	{
		tvec = V0 - RayOrigin;
		dot = -dot;
	}
	if ( dot < 0.0001 )
		return false;

	double u = tvec.dot(pvec);
	// check finite plane.
	if ( u < 0.0 || u > dot )
		return false;

	CVector3 qvec(tvec);
	qvec.cross(edge1);
	double v = RayDirection.dot(qvec);
	// check finite plane.
	if ( v < 0.0 || ( u + v ) > dot ) return false;

	double t = edge2.dot(qvec);
	double flnvDet = 1.0 / dot;
	t *= flnvDet;
	u *= flnvDet;
	v *= flnvDet;

	// calc cross point
	for ( int i = 0; i < 3; ++i )
	{
		cross[i] = RayOrigin[i] + RayDirection[i] * t;
	}

	return true;
}

CVector3::CVector3() : x(0.0), y(0.0), z(0.0)
{
}

CVector3::CVector3(const double& x, const double& y, const double& z)
{
	this->x = x;
	this->y = y;
	this->z = z;
}

CVector3::CVector3(const double coordinate[enAxisNum])
{
	this->x = coordinate[enAxisX];
	this->y = coordinate[enAxisY];
	this->z = coordinate[enAxisZ];
}

CVector3::CVector3(const CVector3& rVector)
{
	this->x = rVector.x;
	this->y = rVector.y;
	this->z = rVector.z;
}

inline const double& CVector3::operator[] (int i) const
{
	if ( i<enAxisX || i>enAxisNum )
	{
	}
	switch ( i )
	{
	case enAxisX: return x;
	case enAxisY: return y;
	case enAxisZ: return z;
	default:
		break;
	}
	return z;
}

inline double& CVector3::operator[] (int i)
{
	if ( i<enAxisX || i>enAxisNum )
	{
	}
	switch ( i )
	{
	case enAxisX: return x;
	case enAxisY: return y;
	case enAxisZ: return z;
	default:
		break;
	}
	return z;
}

inline CVector3& CVector3::operator= (const CVector3& rVector)
{
	x = rVector.x;
	y = rVector.y;
	z = rVector.z;
	return *this;
}

inline bool CVector3::operator== (const CVector3& rVector) const
{
	return (
		fabs(x - rVector.x) < CDgnMathFunc::Get_NormalZero() &&
		fabs(y - rVector.y) < CDgnMathFunc::Get_NormalZero() &&
		fabs(z - rVector.z) < CDgnMathFunc::Get_NormalZero() );
}

inline bool CVector3::operator!= (const CVector3& rVector) const
{
	return (
		fabs(x - rVector.x) > CDgnMathFunc::Get_NormalZero() ||
		fabs(y - rVector.y) > CDgnMathFunc::Get_NormalZero() ||
		fabs(z - rVector.z) > CDgnMathFunc::Get_NormalZero() );
}

inline CVector3 CVector3::operator+ (const CVector3& rVector) const
{
	return CVector3(x + rVector.x, y + rVector.y, z + rVector.z);
}

inline CVector3 CVector3::operator- (const CVector3& rVector) const
{
	return CVector3(x - rVector.x, y - rVector.y, z - rVector.z);
}

inline CVector3 CVector3::operator* (double fScalar) const
{
	return CVector3(fScalar*x, fScalar*y, fScalar*z);
}

inline CVector3 CVector3::operator/ (double fScalar) const
{
	if ( fabs(fScalar) < CDgnMathFunc::Get_NormalZero() );
	return CVector3(x/fScalar, y/fScalar, z/fScalar);
}

inline CVector3 CVector3::operator* (const CVector3& rVector) const
{
	return CVector3(x * rVector.x, y * rVector.y, z * rVector.z);
}

inline CVector3 CVector3::operator/ (const CVector3& rVector) const
{
	if ( fabs(rVector.x) < CDgnMathFunc::Get_NormalZero() );
	if ( fabs(rVector.y) < CDgnMathFunc::Get_NormalZero() );
	if ( fabs(rVector.z) < CDgnMathFunc::Get_NormalZero() );
	return CVector3(x / rVector.x, y / rVector.y, z / rVector.z);
}

inline CVector3 CVector3::operator- () const
{
	return CVector3(-x, -y, -z);
}

inline CVector3& CVector3::operator+= (const CVector3& rVector)
{
	x += rVector.x;
	y += rVector.y;
	z += rVector.z;
	return *this;
}

inline CVector3& CVector3::operator-= (const CVector3& rVector)
{
	x -= rVector.x;
	y -= rVector.y;
	z -= rVector.z;
	return *this;
}

inline CVector3& CVector3::operator*= (double fScalar)
{
	x *= fScalar;
	y *= fScalar;
	z *= fScalar;
	return *this;
}

inline CVector3& CVector3::operator/= (double fScalar)
{
	if ( fabs(fScalar) < CDgnMathFunc::Get_NormalZero() );
	x /= fScalar;
	y /= fScalar;
	z /= fScalar;
	return *this;
}

inline CVector3& CVector3::operator*= (const CVector3& rVector)
{
	x *= rVector.x;
	y *= rVector.y;
	z *= rVector.z;
	return *this;
}

inline CVector3& CVector3::operator/= (const CVector3& rVector)
{
	if ( fabs(rVector.x) < CDgnMathFunc::Get_NormalZero() );
	if ( fabs(rVector.y) < CDgnMathFunc::Get_NormalZero() );
	if ( fabs(rVector.z) < CDgnMathFunc::Get_NormalZero() );
	x /= rVector.x;
	y /= rVector.y;
	z /= rVector.z;
	return *this;
}

void CVector3::GetCoord(double& x, double& y, double& z)
{
	x = this->x;
	y = this->y;
	z = this->z;
}

void CVector3::GetCoord(double coord[3])
{
	coord[0] = this->x;
	coord[1] = this->y;
	coord[2] = this->z;
}

bool CVector3::isZero() const
{
	return squaredLength() == 0.0 ? true : false;
}

bool CVector3::isUnit() const
{
	return squaredLength() == 1.0 ? true : false;
}

double CVector3::length() const
{
	return sqrt(squaredLength());
}

double CVector3::squaredLength() const
{
	return x * x + y * y + z * z;
}

double CVector3::dot(const CVector3& rVector) const
{
	return x * rVector.x + y * rVector.y + z * rVector.z;
}

CVector3 CVector3::cross(const CVector3& rVector) const
{
	return CVector3(y * rVector.z - z * rVector.y,
		z * rVector.x - x * rVector.z,
		x * rVector.y - y * rVector.x);
}

double CVector3::normalize()
{
	double fLength = length();
	if ( fLength > CDgnMathFunc::Get_NormalZero() )
	{
		float fInvLength = 1.0 / fLength;
		x *= fInvLength;
		y *= fInvLength;
		z *= fInvLength;
	}
	else
	{
		fLength = 0.0;
	}

	return fLength;
}

CVector3 CVector3::normalizecross(const CVector3& rVector) const
{
	CVector3 rCross = cross(rVector);
	rCross.normalize();
	return rCross;
}

CVector3 CVector3::minVec(const CVector3 &v) const
{
	return CVector3(__min(v.x, x), __min(v.y, y), __min(v.z, z));
}

CVector3 CVector3::maxVec(const CVector3 &v) const
{
	return CVector3(__min(v.x, x), __min(v.y, y), __min(v.z, z));
}
