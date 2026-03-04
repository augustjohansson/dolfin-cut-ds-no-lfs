#ifndef DOLFIN_SIMPLEX_TOOLS_H
#define DOLFIN_SIMPLEX_TOOLS_H

#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <iomanip>

//#include </home/august/dev/fenics-dev/logg-dolfin/dolfin/dolfin.h>
//#include </home/fenics/shared/dolfin-test/dolfin/dolfin.h>
//#include </home/august/dev/fenics-dev/env_dolfin_logg/logg-dolfin/dolfin/dolfin.h>
//#include </home/fenics/shared/dolfin/dolfin/dolfin.h>
#include "../dolfin.h"

namespace tools
{
#define PPause {char dummycharXohs5su8='a';std::cout<<"\n% Pause: "<<__FILE__<<" line "<<__LINE__<<" function "<<__FUNCTION__<<std::endl;std::cin>>dummycharXohs5su8;}

  using namespace dolfin;

  //-----------------------------------------------------------------------------
  inline std::vector<Point> convert(const double* x,
				    std::size_t tdim,
				    std::size_t gdim)
  {
    std::vector<Point> s(tdim + 1);
    for (std::size_t t = 0; t < tdim + 1; ++t)
      for (std::size_t d = 0; d < gdim; ++d)
	s[t][d] = x[gdim*t + d];
    return s;
  }
  //-----------------------------------------------------------------------------
  inline std::vector<std::vector<Point>> convert(const std::vector<double>& triangulation,
						 std::size_t tdim,
						 std::size_t gdim)
  {
    const std::size_t offset = (tdim + 1)*gdim;
    const std::size_t N = triangulation.size() / offset;
    std::vector<std::vector<Point>> simplices(N);

    for (std::size_t k = 0; k < N; ++k)
    {
      const double* x = triangulation.data() + k*offset;
      simplices[k] = convert(x, tdim, gdim);
    }

    return simplices;
  }

  //-----------------------------------------------------------------------------
  inline std::vector<double> convert(const std::vector<Point>& simplex,
				     std::size_t tdim,
				     std::size_t gdim)
  {
    std::vector<double> x((tdim + 1)*gdim);

    for (std::size_t i = 0; i < tdim + 1; ++i)
      for (std::size_t d = 0; d < gdim; ++d)
	x[i*gdim + d] = simplex[i][d];

    return x;
  }

  //-----------------------------------------------------------------------------
  inline std::vector<double> convert(const std::vector<std::vector<Point>>& simplices,
				     std::size_t tdim,
				     std::size_t gdim)
  {
    const std::size_t offset = (tdim + 1)*gdim;
    const std::size_t N = simplices.size();
    std::vector<double> triangulation(N*offset);

    for (std::size_t i = 0; i < N; ++i)
    {
      const std::vector<double> x = convert(simplices[i], tdim, gdim);
      for (std::size_t j = 0; j < x.size(); ++j)
	triangulation[i*offset + j] = x[j];
    }

    return triangulation;
  }
  //-----------------------------------------------------------------------------
  inline std::vector<Point> convert(const Cell& cell)
  {
    const std::size_t tdim = cell.mesh().topology().dim();
    std::vector<Point> simplex(tdim + 1);
    const MeshGeometry& geometry = cell.mesh().geometry();
    const unsigned int* vertices = cell.entities(0);
    for (std::size_t j = 0; j < tdim + 1; ++j)
      simplex[j] = geometry.point(vertices[j]);
    return simplex;
  }

  //-----------------------------------------------------------------------------
  // typedef std::vector<Point> Simplex;
  // typedef std::vector<Simplex> Polyhedron;
  inline bool tdimcheck(const std::vector<std::vector<dolfin::Point>>& polygon)
  {
    if (polygon.size() == 0) return false;

    const std::size_t tdimtmp = polygon[0].size();
    for (std::size_t i = 1; i < polygon.size(); ++i)
      if (polygon.at(i).size() != tdimtmp)
	return false;
    return true;
  }

  //-----------------------------------------------------------------------------
  inline bool tdimcheck(const std::vector<std::vector<std::vector<dolfin::Point>>>& pvec)
  {
    if (pvec.size() == 0) return false;
    for (std::size_t i = 0; i < pvec.size(); ++i)
      if (pvec[i].size() == 0)
	return false;

    const std::size_t tdimtmp = pvec[0][0].size();
    for (std::size_t i = 1; i < pvec.size(); ++i)
      for (std::size_t j = 0; j < pvec[i].size(); ++j)
	if (pvec.at(i).at(j).size() != tdimtmp)
	  return false;
    return true;
  }


  //-----------------------------------------------------------------------------
  // display quadrature_rule
  // recall typedef std::pair<std::vector<double>, std::vector<double> > quadrature_rule;
  inline void cout_qr(const std::pair<std::vector<double>, std::vector<double> >& qr,
		      std::string color="'b'",
		      std::size_t markersize=16)
  {
    std::cout << "% "<<__FUNCTION__<<' '<<qr.second.size()<<std::endl;
    double vol = 0.0;
    for (std::size_t i = 0; i < qr.second.size(); ++i)
    {
      vol += qr.second[i];
      std::stringstream ss;
      if (qr.second[i] > 0)
	ss<<"'color',"<<color<<",'marker','.','markersize',"<<markersize;
      else
	ss<<"'color',"<<color<<",'marker','o','markersize',"<<markersize-10;
      //std::cout << "plot("<<qr.first[2*i]<<','<<qr.first[2*i+1]<<','<<marker<<"); % "<<qr.second[i]<<' '<<i<<std::endl;
      // if (qr.first[2*i+1]<qr.first[2*i])
      std::cout << "plot("<<qr.first[2*i]<<','<<qr.first[2*i+1]<<','<<ss.str()<<"); % "<<std::setprecision(std::numeric_limits<long double>::digits10+2)<<qr.second[i]<<' '<<i<<std::endl;
    }
    std::cout<<"% net vol " << vol<<std::endl;
  }
  
  //-----------------------------------------------------------------------------
  inline void cout_qr(const std::vector<std::pair<std::vector<double>, std::vector<double>>>& qr,
                      const std::vector<std::string>& color={},
                      std::size_t markersize=16)
  {
    auto c = color;
    if (color.size() == 0)
      for (std::size_t i = 0; i < qr.size(); ++i)
        c.push_back("'b'");
      
    for (std::size_t i = 0; i < qr.size(); ++i)
      cout_qr(qr[i], c[i], markersize);
  }
  
  //-----------------------------------------------------------------------------
  inline void cout_normals(const std::vector<double>& n)
  {
    for (std::size_t i = 0; i < n.size()/2; ++i)
      std::cout << i << ":  "<<n[2*i]<<' '<<n[2*i+1]<<std::endl;
  }

  //-----------------------------------------------------------------------------
  inline double area(const std::pair<std::vector<double>, std::vector<double> >& qr)
  {
    double a = 0;
    for (std::size_t i = 0; i < qr.second.size(); ++i)
      a += qr.second[i];
    return a;
  }
  //-----------------------------------------------------------------------------
  inline std::string weights(const std::pair<std::vector<double>, std::vector<double> >& qr)
  {
    std::stringstream ss;
    for (std::size_t i = 0; i < qr.second.size(); ++i)
      ss << qr.second[i]<<' ';
    return ss.str();
  }
  //-----------------------------------------------------------------------------
  inline double volume_of_cut_part(const MultiMesh& multimesh,
				   const std::size_t p)
  {
    double volume = 0.0;

    // Sum volume of uncut cells (from cell.volume)
    {
      const auto& cells = multimesh.uncut_cells(p);
      for (auto it = cells.begin(); it != cells.end(); ++it)
      {
	const Cell cell(*multimesh.part(p), *it);
	volume += cell.volume();
      }
    }

    // Sum volume of cut cells (from quadrature rules)
    {
      const auto& cells = multimesh.cut_cells(p);
      for (auto it = cells.begin(); it != cells.end(); ++it)
      {
	const auto& qr = multimesh.quadrature_rules_cut_cells(p, *it);
	for (std::size_t i = 0; i < qr.second.size(); ++i)
	  volume += qr.second[i];
      }
    }

    return volume;
  }
  //-----------------------------------------------------------------------------
  inline double area_of_cut_part(const MultiMesh& multimesh,
				 const std::size_t cut_part)
  {
    double a = 0.0;
    const auto& quadrature_rules = multimesh.quadrature_rules_interface(cut_part);

    // Get collision map
    const auto& cmap = multimesh.collision_map_cut_cells(cut_part);
    for (auto it = cmap.begin(); it != cmap.end(); ++it)
    {
      const unsigned int cut_cell_index = it->first;
      const auto& cutting_cells = it->second;

      // Iterate over cutting cells
      for (auto jt = cutting_cells.begin(); jt != cutting_cells.end(); jt++)
      {
	// Get quadrature rule for interface part defined by
	// intersection of the cut and cutting cells
	const std::size_t k = jt - cutting_cells.begin();
	dolfin_assert(k < quadrature_rules.at(cut_cell_index).size());
	const auto& qr = quadrature_rules.at(cut_cell_index)[k];
	a += area(qr);
      }
    }

    return a;
  }

  //-----------------------------------------------------------------------------
  // this sorts such that a >= b >= c
  template<class T>
  inline void sort3(T &a, T &b, T &c)
  {
    if (b>a) std::swap(b,a);
    if (c>b) std::swap(c,b);
    if (b>a) std::swap(b,a);
  }

  //-----------------------------------------------------------------------------
  inline double Heron(double a, double b, double c)
  {
    sort3(a,b,c);
    const double s2 = (a+(b+c))*(c-(a-b))*(c+(a-b))*(a+(b-c));
    if (s2 < 0)
    {
      std::cout << "Heron error, negative sqrt: " << s2 << " is to be replaced with 0" << std::endl;
      if (std::abs(s2) < DOLFIN_EPS)
  	return 0;
      else
  	exit(1);
    }
    return 0.25*std::sqrt( (a+(b+c))*(c-(a-b))*(c+(a-b))*(a+(b-c)) );
  }

  //-----------------------------------------------------------------------------
  inline double area(const std::vector<dolfin::Point> &simplex)
  {
    if (simplex.size() == 3)
    {
      const dolfin::Point et=simplex[1]-simplex[0];
      const dolfin::Point es=simplex[2]-simplex[0];
      return Heron(et.norm(), es.norm(), (et-es).norm());
      // double a[] = {simplex[0].x(), simplex[0].y()};
      // double b[] = {simplex[1].x(), simplex[1].y()};
      // double c[] = {simplex[2].x(), simplex[2].y()};
      // return 0.5*std::abs(orient2d(a,b,c));
    }
    else if (simplex.size() == 2)
    {
      return (simplex[0]-simplex[1]).norm();
    }
    else if (simplex.size() == 1)
    {
      return 0.;
    }
    else {
      std::cout << "error simplex size = " << simplex.size();
      PPause;
      return -9e99;
    }

    // if (simplex.size() == 3)
    // {
    //   double a[2]={simplex[0][0],simplex[0][1]};
    //   double b[2]={simplex[1][0],simplex[1][1]};
    //   double c[2]={simplex[2][0],simplex[2][1]};
    //   return 0.5*orient2d(a,b,c);
    // }
    // else if (simplex.size() == 2)
    // {
    //   return (simplex[0]-simplex[1]).norm();
    // }
    // else
    // {
    //   PPause;
    //   return -9e99;
    // }
  }

  //-----------------------------------------------------------------------------
  inline std::string drawtriangle(const std::vector<dolfin::Point> &simplex,
				  const std::string& color = "[0.6 0.7 0.8]",
				  bool matlab=true,
				  bool threeD=true)
  {
    std::stringstream ss; ss.precision(std::numeric_limits<long double>::digits10+2);
    if (simplex.size() == 4)
    {
      if (matlab)
	ss << "drawtet(";
      else { PPause; }
      ss<< "["<<simplex[0][0]<<','<<simplex[0][1]<<','<<simplex[0][2]<<"],"
	<< "["<<simplex[1][0]<<','<<simplex[1][1]<<','<<simplex[1][2]<<"],"
	<< "["<<simplex[2][0]<<','<<simplex[2][1]<<','<<simplex[2][2]<<"],"
	<< "["<<simplex[3][0]<<','<<simplex[3][1]<<','<<simplex[3][2]<<"]";
      if (matlab)
	ss << ","<<color<<");";
    }
    else if (simplex.size() == 3)
    {
      if (matlab)
	ss << "drawtriangle(";
      else
	ss << "drawtriangle2(";
      if (threeD)
      {
	ss<< "["<<simplex[0][0]<<','<<simplex[0][1]<<','<<simplex[0][2]<<"],"
	  << "["<<simplex[1][0]<<','<<simplex[1][1]<<','<<simplex[1][2]<<"],"
	  << "["<<simplex[2][0]<<','<<simplex[2][1]<<','<<simplex[2][2]<<"]";
      }
      else
      {
      ss<< "["<<simplex[0][0]<<','<<simplex[0][1]<<"],"
	<< "["<<simplex[1][0]<<','<<simplex[1][1]<<"],"
	<< "["<<simplex[2][0]<<','<<simplex[2][1]<<"]";
      }
      if (matlab)
	ss << ","<<color<<");";
      else {
	ss << ",color=" << color << ',' //<< ");";
	   << "plt=plt,axis=gca()"
	   << ");";
      }
    }
    else if (simplex.size() == 2)
    {
      if (threeD)
      {
	ss << "drawline([" << simplex[0][0] << ',' << simplex[0][1]<<','<<simplex[0][2] << "],"
	   <<  "[" << simplex[1][0] << ',' << simplex[1][1]<<','<<simplex[1][2] << "],";
      }
      else {
	ss << "drawline([" << simplex[0][0] << ',' << simplex[0][1] << "],"
	   <<  "[" << simplex[1][0] << ',' << simplex[1][1] << "],";
      }
      if (matlab)
	ss << color<<",1,1,15);";
      else
	ss << "plt=plt,color="<< color << ",linewidth=5.0);";
    }
    else if (simplex.size() == 1)
    {
      ss << "plot("<<simplex[0][0]<<','<<simplex[0][1]<<',';
      if (matlab)
	ss << "'k.','markersize',15);";
      else {
	PPause;
	dolfin_assert(false); // /not implemented
      }
    }
    else {
      std::cout << "simplex size to plot is " << simplex.size() << std::endl;
      PPause;
      dolfin_assert(false);
    }
    return ss.str();
  }

  //-----------------------------------------------------------------------------
  inline std::string drawtriangle3(const std::vector<dolfin::Point> &simplex,
				   const std::string& color = "[0.6 0.7 0.8]")
  {
    return drawtriangle(simplex, color, true, true);
  }

  //-----------------------------------------------------------------------------
  inline std::string drawtriangle(const dolfin::Cell &cell,
				  const std::string& color = "[0.6 0.7 0.8]",
				  bool matlab=true)
  {
    const std::vector<Point> s = convert(cell);
    return drawtriangle(s, color, matlab);
  }

  //-----------------------------------------------------------------------------
  inline std::string drawtriangle(const std::vector<double>& s,
				  const std::string& color = "[0.6 0.7 0.8]",
				  bool matlab = false)
  {
    std::vector<dolfin::Point> pp(s.size() / 2);
    for (std::size_t i = 0; i < pp.size(); ++i)
    {
      pp[i][0] = s[2*i];
      pp[i][1] = s[2*i+1];
    }
    return drawtriangle(pp, color, matlab);

    // std::vector<dolfin::Point> ss(3);
    // ss[0] = dolfin::Point(s[0],s[1]);
    // ss[1] = dolfin::Point(s[2],s[3]);
    // ss[2] = dolfin::Point(s[4],s[5]);
    // return drawtriangle(ss, color);
  }

  //-----------------------------------------------------------------------------
  inline std::string drawtriangle(const dolfin::Point a,
				  const dolfin::Point b,
				  const dolfin::Point c,
				  const std::string color = "[0.6 0.7 0.8]",
				  bool matlab = false)
  {
    std::vector<dolfin::Point> t = {{ a, b, c}};
    return drawtriangle(t, color);
  }

  //------------------------------------------------------------------------------
  inline std::string drawline(const dolfin::Point a,
			      const dolfin::Point b,
			      const std::string color = "'b'",
			      bool matlab = true)
  {
    return drawtriangle({{a,b}},color,matlab);
  }

  //-----------------------------------------------------------------------------
  inline std::string drawtriangle(const std::vector<std::vector<dolfin::Point> >& polyhedron,
				  const std::string color = "[0.6 0.7 0.8]")
  {
    std::stringstream ss;
    for (const std::vector<dolfin::Point>& simplex : polyhedron)
      ss << drawtriangle(simplex, color);
    return ss.str();
  }

  //-----------------------------------------------------------------------------
  inline std::string matlabplot2(const dolfin::Point p,
				 const std::string m="'k.','markersize',14")
  {
    std::stringstream ss; ss.precision(std::numeric_limits<long double>::digits10+2);
    ss<<"plot3("<<p[0]<<','<<p[1]<<','<<m<<");";
    return ss.str();
  }
  inline std::string matlabplot3(const dolfin::Point p,
  				const std::string m="'k.','markersize',14")
  {
    std::stringstream ss; ss.precision(std::numeric_limits<long double>::digits10+2);
    ss<<"plot3("<<p[0]<<','<<p[1]<<','<<p[2]<<','<<m<<");";
    return ss.str();
  }

  //-----------------------------------------------------------------------------
  inline std::string plot(const dolfin::Point p,
			  const std::string m="'k.','markersize',14")
  {
    return matlabplot3(p,m);
  }

  inline std::string plot2(const dolfin::Point p,
			  const std::string m="'k.','markersize',14")
  {
    return matlabplot2(p,m);
  }

  inline std::string plot3(const dolfin::Point p,
			  const std::string m="'k.','markersize',14")
  {
    return matlabplot3(p,m);
  }


  //-----------------------------------------------------------------------------
  inline std::string drawarrow(const dolfin::Point v1,
			       const dolfin::Point v2,
			       const std::string& color = "'b'")
  {
    std::stringstream ss;ss.precision(std::numeric_limits<long double>::digits10+2);
    ss << "drawarrow([" << v1[0] << ' '<<v1[1] <<"],[" << v2[0]<<' '<<v2[1] << "], "<< color << ");";
    return ss.str();
    // const dolfin::Point v = v2-v1;
    // const Point ones(0,0,1);
    // Point n = ones.cross(v);
    // if (n.norm() < 1e-5) {
    //   const Point ones(0,1,0);
    //   n = ones.cross(v);
    // }
    // const double a = 0.03*norm(v);
    // n /= n.norm();
    // drawline(v1, v2);
    // drawline(v2, v1 + 0.8 * v + a * n);
    // drawline(v2, v1 + 0.8 * v - a * n);
  }

  //void Pause() { char apa; std::cin >> apa; }
  //-----------------------------------------------------------------------------
  inline std::string drawcircle(const Point& p,
				double r,
				const std::string& color= "'b'")
  {
    std::stringstream ss;
    ss.precision(std::numeric_limits<long double>::digits10+2);
    ss << "drawcircle([" << p[0]<<","<<p[1] <<"]," << r << "," << color << ");";
    return ss.str();
  }

  //-----------------------------------------------------------------------------
  template<class U=std::size_t, class T=double>
  inline void dolfin_write_medit(const std::string &filename,
				 const dolfin::Mesh& mesh,
				 const U t = 0,
				 const std::vector<T>* u = 0)
  {
    const std::size_t tdim = mesh.topology().dim();
    const std::size_t gdim = mesh.geometry().dim();

    std::stringstream ss;
    ss<<filename<<"."<<t<<".mesh";
    std::ofstream file(ss.str().c_str());
    if (!file.good()) { std::cout << "sth wrong with the file " << ss.str()<<'\n'; exit(0); }
    file.precision(13);
    // write vertices
    const std::size_t nno = mesh.num_vertices();
    file << "MeshVersionFormatted 1\nDimension\n"<<gdim<<"\nVertices\n"
	 << nno<<'\n';
    const std::vector<double>& coords = mesh.coordinates();
    for (std::size_t i = 0; i < nno; ++i) {
      for (std::size_t d = 0; d < gdim; ++d)
	file << coords[gdim*i+d]<<' ';
      file <<" 1\n";
    }
    // write connectivity
    const std::size_t nel = mesh.num_cells();
    file << ((tdim == 2) ? "Triangles\n" : "Tetrahedra\n")
	 << nel <<'\n';
    const std::vector<unsigned int>& cells = mesh.cells();
    for (std::size_t e = 0; e < nel; ++e) {
      for (std::size_t t = 0; t < tdim+1; ++t)
	file << cells[(tdim+1)*e+t]+1<<' ';
      if (u and u->size()==nel)
	file << " " << (*u)[e] << '\n';
      else
	file << " 1\n";
    }
    file.close();

    if (u)
    {
      std::stringstream ss;
      ss<<filename<<"."<<t<<".bb";
      std::ofstream file(ss.str().c_str());
      if (!file.good()) { std::cout << "sth wrong with the file " << ss.str()<<'\n'; exit(0); }
      file.precision(13);
      const std::size_t nno = mesh.num_vertices();
      const std::size_t nel = mesh.num_cells();
      // Write data (node or element based).
      if (u->size()==nno)
	file << "3 1 " << nno << " 2\n";
      else if (u->size()==nel)
	file << "3 1 " << nel << " 1\n";
      else
      {
	std::cout<<"\n\nstrange sizes u="<<u->size() << " nno=" << nno<<" nel="<<nel << ". Writing available data anyway\n"<<std::endl;
	file << "3 1 " << u->size() << " 1\n";
      }

      // Writing
      for (std::size_t i = 0; i < u->size(); ++i)
	file << (*u)[i] << '\n';
      file.close();
    }
  }
  //-----------------------------------------------------------------------------
  inline void writematlab(const std::string& filename,
			  const dolfin::MultiMesh& mm,
			  const bool flat=true,
			  const int t=0)
  {
    const bool threeD=!flat;
    std::stringstream ss;
    ss<<filename<<"_"<<t<<".m";
    std::ofstream file(ss.str().c_str());
    if (!file.good()) { std::cout << "sth wrong with the file " << ss.str()<<'\n'; exit(0); }
    file.precision(13);

    for (std::size_t i=0; i<mm.num_parts(); ++i)
    {
      for (std::size_t e=0; e<mm.part(i)->num_cells(); ++e)
      {
	const Cell cell(*mm.part(i), e);
	std::vector<Point> s = convert(cell);
	for (Point p: s)
	  p[2] = (flat) ? 0 : i;
	file << drawtriangle(s,"[0.6 0.7 0.8]",true,threeD);
      }
    }
    if (flat)
      file << ";axis square; axis tight; xlabel x; ylabel y; view(2);\n";
    else
      file << ";axis square; axis tight; xlabel x; ylabel y; view(-28,10);\n";
    file.close();
  }

  //-----------------------------------------------------------------------------
  inline void dolfin_write_medit_parts(const std::string& filename,
				       const dolfin::MultiMesh& mm)
  {
    for (std::size_t p=0; p<mm.num_parts(); ++p)
    {
      const auto uncut = mm.uncut_cells(p);
      const auto cut = mm.cut_cells(p);
      const auto covered = mm.covered_cells(p);
      std::vector<int> marker(mm.part(p)->num_cells(),-1);
      // 0: uncut   = cell not colliding with any higher domain
      // 1: cut     = cell colliding with some higher boundary and is not covered
      // 2: covered = cell colliding with some higher domain but not its boundary
      for (const auto c: uncut) marker[c] = 0;
      for (const auto c: cut) marker[c] = 1;
      for (const auto c: covered) marker[c] = 2;

      dolfin_write_medit(filename,*mm.part(p),p,&marker);
    }
  }

  //-----------------------------------------------------------------------------
  inline void dolfin_write_medit(const std::string& filename,
				 const dolfin::MultiMesh& mm,
				 //const std::vector<std::vector<double>> *u=0,
				 const int t=0)
  {
    std::stringstream ss;
    ss<<filename<<"."<<t<<".mesh";
    std::ofstream file(ss.str().c_str());
    if (!file.good()) { std::cout << "sth wrong with the file " << ss.str()<<'\n'; exit(0); }
    file.precision(13);

    assert(mm.num_parts()>0);
    const std::size_t tdim = mm.part(0)->topology().dim();
    const std::size_t gdim = mm.part(0)->geometry().dim();

    // write vertices
    std::size_t nno=0;
    for (std::size_t i=0; i<mm.num_parts(); ++i)
      nno += mm.part(i)->num_vertices();
    file << "MeshVersionFormatted 1\nDimension\n"<<gdim<<"\nVertices\n"
	 << nno<<'\n';
    for (std::size_t i=0; i<mm.num_parts(); ++i) {
      const std::vector<double>& coords = mm.part(i)->coordinates();
      for (std::size_t j=0; j<mm.part(i)->num_vertices(); ++j) {
	for (std::size_t d=0; d<gdim; ++d)
	  file << coords[gdim*j+d]<<' ';
	file << i+1 << '\n';
      }
    //   {
    // 	if (gdim == 3)
    // 	  file << coords[2*j]<<' '<<coords[2*j+1]<<' '<<(double)i/mm.num_parts()<<' '<<i+1<<'\n';
    // 	else if (gdim == 2)
    // 	file << coords[2*j]<<' '<<coords[2*j+1]<<' '<<i+1<<'\n';
    // 	else { PPause; }
    //   }
    }

    // write connectivity
    std::size_t nel=0;
    for (std::size_t i=0; i<mm.num_parts(); ++i)
      nel += mm.part(i)->num_cells();
    file << ((tdim == 2) ? "Triangles\n" : "Tetrahedra\n")
	 << nel <<'\n';
    const std::size_t nne=tdim+1;
    std::size_t offset = 0;
    for (std::size_t i=0; i<mm.num_parts(); ++i) {
      const std::vector<unsigned int>& cells = mm.part(i)->cells();
      offset += (i == 0) ? 0 : mm.part(i-1)->num_vertices();
      for (std::size_t e = 0; e < mm.part(i)->num_cells(); ++e) {
	for (std::size_t j=0; j<nne; ++j)
	  file << cells[nne*e+j]+offset+1<<' ';
	file << i+1<<'\n';
      }

      //file << cells[3*e]+offset+1<<' '<<cells[3*e+1]+offset+1<<' '<<cells[3*e+2]+offset+1<<' '<<i+1<<'\n';
      // for (const auto e: mm.uncut_cells(i))
      // 	file << cells[3*e]+offset+1<<' '<<cells[3*e+1]+offset+1<<' '<<cells[3*e+2]+offset+1<<' '<<i+1<<'\n';
      // for (const auto e: mm.cut_cells(i))
      // 	file << cells[3*e]+offset+1<<' '<<cells[3*e+1]+offset+1<<' '<<cells[3*e+2]+offset+1<<' '<<i+1<<'\n';

    }
    file.close();

    {
      std::stringstream ss;
      ss<<filename<<"."<<t<<".bb";
      std::ofstream file(ss.str().c_str());
      if (!file.good()) { std::cout << "sth wrong with the file " << ss.str()<<'\n'; exit(0); }
      file.precision(13);
      file << "3 1 " << nel << " 1\n";

      for (std::size_t i=0; i<mm.num_parts(); ++i)
	for (std::size_t j=0; j<mm.part(i)->num_cells(); ++j)
	  file << i+1 <<'\n';
      file.close();
    }

  }


  //-----------------------------------------------------------------------------
  // Hack to write vtu file
  inline void write_vtu_hack(const std::string& filename,
			     const dolfin::Mesh& mesh,
			     std::size_t i = 0)
  {
    std::stringstream ss;
    ss << i;
    std::ofstream fp(filename+ss.str()+".vtu");

    const std::size_t num_vertices = mesh.num_vertices();
    const std::size_t num_cells = mesh.num_cells();

    fp << "<?xml version=\"1.0\"?>\n"
       << "<VTKFile type=\"UnstructuredGrid\"  version=\"0.1\"  >\n"
       << "<UnstructuredGrid>\n"
       << "<Piece  NumberOfPoints=\"" << num_vertices << "\" NumberOfCells=\"" << num_cells << "\">\n"
       << "<Points>\n"
       << "<DataArray  type=\"Float64\"  NumberOfComponents=\"3\"  format=\"ascii\">";

    // vertices
    for (dolfin::VertexIterator vertex(mesh); !vertex.end(); ++vertex) {
      for (int d = 0; d < 2; ++d) // dimension
	fp << vertex->x(d) << ' ';
      fp << "0   "; // always write 3d
    }
    fp << "</DataArray>\n"
       << "</Points>\n";

    // cells
    fp << "<Cells>\n"
       << "<DataArray  type=\"UInt32\"  Name=\"connectivity\"  format=\"ascii\">";
    const std::vector<unsigned int>& cells = mesh.cells();
    for (std::size_t e = 0; e < num_cells; ++e) {
      // tets:
      //fp << cells[4*e] << ' ' << cells[4*e+1] << ' ' << cells[4*e+2] << ' ' << cells[4*e+3] << "  ";
      // tris:
      fp << cells[3*e]<<' '<<cells[3*e+1]<<' '<<cells[3*e+2]<<"  ";
    }
    fp << "</DataArray>\n";

    // offset
    fp << "<DataArray  type=\"UInt32\"  Name=\"offsets\"  format=\"ascii\">";
    for (std::size_t e = 0, offset=3; e < num_cells; ++e, offset += 3) // offset is 3 or 4
      fp << offset << ' ';
    fp << "</DataArray>\n";

    // types
    const std::size_t vtk_element_type = 5; // tet=10, tri=5
    fp << "<DataArray  type=\"UInt8\"  Name=\"types\"  format=\"ascii\">";
    for (std::size_t e = 0; e < num_cells; ++e)
      fp << vtk_element_type << ' ';
    fp << "</DataArray>\n"
       << "</Cells>\n";

    // data
    fp.precision(16);
    const std::size_t size = num_vertices;
    std::vector<double> values(size, i);
    //u.compute_vertex_values(values, mesh);
    const std::string encode_string = "ascii";

    // // write velocity
    // const std::string velocity_name = u.name() + "_velocity";
    // fp << "<PointData>\n"
    //    << "<DataArray  type=\"Float64\"  Name=\"" << velocity_name << "\"  NumberOfComponents=\"3\" format=\""<< encode_string <<"\">";
    // for (VertexIterator vertex(mesh); !vertex.end(); ++vertex)
    // {
    //   for (std::size_t i = 0; i < 3; ++i) // Only write 3 components!
    // 	fp << values[vertex->index() + i*num_vertices] << " ";
    //   fp << " ";
    // }
    // fp << "</DataArray>\n";

    // // write pressure
    // const std::string pressure_name = u.name() + "_pressure";
    // fp << "<DataArray  type=\"Float64\"  Name=\"" << pressure_name << "\"  NumberOfComponents=\"1\" format=\""<< encode_string <<"\">";
    // for (VertexIterator vertex(mesh); !vertex.end(); ++vertex)
    //   fp << values[vertex->index() + 3*num_vertices] << ' ';
    // fp << "</DataArray>\n"
    //    << "</PointData>\n";

    const std::string name = "data_part_"+ss.str();
    fp << "<PointData>\n"
      //<< "<DataArray  type=\"Float64\"  Name=\"" << name << "\"  NumberOfComponents=\"1\" format=\""<< encode_string <<"\">";
       << "<DataArray  type=\"Float64\"  Name=\"" << name << "\" format=\""<< encode_string <<"\">";
    for (dolfin::VertexIterator vertex(mesh); !vertex.end(); ++vertex)
      fp << values[vertex->index()] << ' ';
    fp << "</DataArray>\n"
       << "</PointData>\n";


    fp << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
    fp.close();
  }

  //-----------------------------------------------------------------------------
  inline std::string zoom(bool matlab=true)
  {
    if (matlab) return "axis equal;";
    else
      return "plt.autoscale(enable=True,axis='both',tight=None);";
  }

  //-----------------------------------------------------------------------------
  inline void writemarkers(const std::string& filename,
			   const MultiMesh& mm,
			   std::size_t step = 0)
  {
    for (std::size_t part = 0; part < mm.num_parts(); ++part)
    {
      std::stringstream ss; ss << part;
      const std::size_t n = mm.part(part)->num_cells();
      std::vector<int> uncut(n, -1), cut(n, -1), covered(n, -1);
      for (const auto c: mm.uncut_cells(part)) uncut[c] = 0;
      for (const auto c: mm.cut_cells(part)) cut[c] = 1;
      for (const auto c: mm.covered_cells(part)) covered[c] = 2;
      dolfin_write_medit(filename+"_uncut"+ss.str(),*mm.part(part),step,&uncut);
      dolfin_write_medit(filename+"_cut"+ss.str(),*mm.part(part),step,&cut);
      dolfin_write_medit(filename+"_covered"+ss.str(),*mm.part(part),step,&covered);
    }
    dolfin_write_medit(filename,mm,step);

  }

  //------------------------------------------------------------------------------
  inline double compute_volume(const MultiMesh& multimesh)
  {
    std::cout << std::endl << __FUNCTION__<< std::endl;

    double volume = 0;
    std::vector<double> all_volumes;

    std::ofstream file("quadrature_volume.txt");
    if (!file.good()) { std::cout << "file not good"<<std::endl; exit(0); }
    file.precision(std::numeric_limits<long double>::digits10+2);

    // Sum contribution from all parts
    std::cout << "Sum contributions"<<std::endl;
    for (std::size_t part = 0; part < multimesh.num_parts(); part++)
    {
      std::cout << "% part " << part;
      double part_volume = 0;
      std::vector<double> status(multimesh.part(part)->num_cells(), 0);

      // Uncut cell volume given by function volume
      double uncut_volume = 0;
      const auto uncut_cells = multimesh.uncut_cells(part);
      for (auto it = uncut_cells.begin(); it != uncut_cells.end(); ++it)
      {
	const Cell cell(*multimesh.part(part), *it);
	file << cell.midpoint().x() << ' ' << cell.midpoint().y()<<' '<<cell.volume() << std::endl;
	volume += cell.volume();
	part_volume += cell.volume();
	uncut_volume += cell.volume();
	status[*it] = 1;
      }

      std::cout << "\t uncut volume "<< uncut_volume<<' ';

      // Cut cell volume given by quadrature rule
      double cut_volume = 0;
      const auto& cut_cells = multimesh.cut_cells(part);
      for (auto it = cut_cells.begin(); it != cut_cells.end(); ++it)
      {
	const auto& qr = multimesh.quadrature_rules_cut_cells(part,*it); PPause; //multimesh.quadrature_rule_cut_cell(part, *it);
	for (std::size_t i = 0; i < qr.second.size(); ++i)
	{
	  file << qr.first[2*i]<<' '<<qr.first[2*i+1]<<' '<<qr.second[i]<<std::endl;
	  volume += qr.second[i];
	  part_volume += qr.second[i];
	  cut_volume += qr.second[i];
	}
	status[*it] = 2;
      }

      std::cout << "\tcut volume " << cut_volume << "\ttotal volume " << part_volume << std::endl;

      all_volumes.push_back(part_volume);

      dolfin_write_medit("status",*multimesh.part(part),part,&status);
    }
    file.close();

    return volume;
  }

  //-----------------------------------------------------------------------------
  inline double compute_volume_overlap(const MultiMesh& multimesh)
  {
    std::cout << std::endl << __FUNCTION__ << std::endl;

    // Mimic MultiMeshAssembler::_assemble_overlap
    double vol = 0;

    // Iterate over parts
    for (std::size_t part = 0; part < multimesh.num_parts(); part++)
    {
      double vol_part = 0;

      // Get quadrature rules
      const auto& quadrature_rules = multimesh.quadrature_rules_overlap(part);PPause;

      // Get collision map
      const auto& cmap = multimesh.collision_map_cut_cells(part);
      // Iterate over all cut cells in collision map
      for (auto it = cmap.begin(); it != cmap.end(); ++it)
      {
	// Get cut cell
	const unsigned int cut_cell_index = it->first;
	const Cell cut_cell(*multimesh.part(part), cut_cell_index);

	// Iterate over cutting cells
	const auto& cutting_cells = it->second;
	for (auto jt = cutting_cells.begin(); jt != cutting_cells.end(); jt++)
	{
	  // Get cutting part and cutting cell
	  const std::size_t cutting_part = jt->first;
	  const std::size_t cutting_cell_index = jt->second;
	  const Cell cutting_cell(*multimesh.part(cutting_part), cutting_cell_index);

	  // Get quadrature rule for interface part defined by
	  // intersection of the cut and cutting cells
	  const std::size_t k = jt - cutting_cells.begin();
	  dolfin_assert(k < quadrature_rules.at(cut_cell_index).size());
	  const auto& qr = quadrature_rules.at(cut_cell_index)[k];

	  // Skip if there are no quadrature points
	  const std::size_t num_quadrature_points = qr.second.size();

	  if (num_quadrature_points > 0)
	  {
	    for (std::size_t i = 0; i < num_quadrature_points; ++i)
	    {
	      vol_part += qr.second[i];
	      vol += qr.second[i];
	    }
	  }
	}
      }
      std::cout << " part " << part << " overlap volume = " << vol_part << std::endl;
    }
    std::cout << " total overlap volume = " << vol << std::endl;
    return vol;
  }

  //------------------------------------------------------------------------------
  inline double compute_interface_area(const MultiMesh& multimesh)
  {
    std::cout << std::endl << __FUNCTION__ << std::endl;

    double area = 0;
    std::vector<double> all_areas;

    std::ofstream file("quadrature_interface.txt");
    if (!file.good()) { std::cout << "file not good"<<std::endl; exit(0); }
    file.precision(std::numeric_limits<long double>::digits10+2);

    // Sum contribution from all parts
    std::cout << "Sum contributions"<<std::endl;
    for (std::size_t part = 0; part < multimesh.num_parts(); part++)
    {
      std::cout << "% part " << part << ' ';
      double part_area = 0;
      const auto& quadrature_rules = multimesh.quadrature_rules_interface(part);PPause;

      // // Uncut cell area given by function area
      // const auto uncut_cells = multimesh.uncut_cells(part);
      // for (auto it = uncut_cells.begin(); it != uncut_cells.end(); ++it)
      // {
      //   const Cell cell(*multimesh.part(part), *it);
      //   area += cell.area();
      // 	//std::cout << std::setprecision(std::numeric_limits<long double>::digits10+2) << cell.area() <<std::endl;
      //   part_area += cell.area();
      // 	status[*it] = 1;
      // 	//file << "0 0 "<< cell.area() << std::endl;
      // }

      // std::cout << "\t uncut area "<< part_area << ' ';


      // Get collision map
      const auto& cmap = multimesh.collision_map_cut_cells(part);
      for (auto it = cmap.begin(); it != cmap.end(); ++it)
      {
	const unsigned int cut_cell_index = it->first;
	const auto& cutting_cells = it->second;

	// Iterate over cutting cells
	for (auto jt = cutting_cells.begin(); jt != cutting_cells.end(); jt++)
	{
	  // Get quadrature rule for interface part defined by
	  // intersection of the cut and cutting cells
	  const std::size_t k = jt - cutting_cells.begin();
	  // std::cout << cut_cell_index << ' ' << k <<' ' << std::flush
	  // 	    << quadrature_rules.size() << ' '
	  // 	    << quadrature_rules.at(cut_cell_index).size() << "   " << std::flush;
	  dolfin_assert(k < quadrature_rules.at(cut_cell_index).size());
	  const auto& qr = quadrature_rules.at(cut_cell_index)[k];
	  std::stringstream ss;
	  for (std::size_t i = 0; i < qr.second.size(); ++i)
	  {
	    file << qr.first[2*i]<<' '<<qr.first[2*i+1]<<' '<<qr.second[i]<<std::endl;
	    //std::cout << qr.second[i]<<' ';
	    area += qr.second[i];
	    part_area += qr.second[i];
	    //std::cout << qr.first[2*i]<<' '<<qr.first[2*i+1]<<std::endl;
	  }
	  //std::cout << std::endl;
	}
      }
      std::cout << "total area " << part_area << std::endl;
      all_areas.push_back(part_area);
    }
    file.close();

    return area;
  }

  //------------------------------------------------------------------------------
  // inline void plot_normals_exterior(const MultiMesh& multimesh)
  // {
  //   std::cout << '%' << __FUNCTION__ << std::endl
  //             << "clf; hold on;\n";
  //   const std::vector<std::string> colors = {{ "'b'", "'g'", "'r'" }};
  //   const std::vector<std::string> marker = {{ "'.'", "'o'", "'x'" }};

  //   for (std::size_t part = 0;
  //        part < std::min<std::size_t>(3, multimesh.num_parts());
  //        part++)
  //     // const std::size_t part = 1;
  //   {
  //     std::cout << "% part " << part << ' ' <<std::endl;
  //     const auto& qr_exterior = multimesh.quadrature_rules_exterior_cut_facets(part);
  //     const auto& facet_normals = multimesh.facet_normals_exterior_cut_facets(part);
  //     for (auto it = qr_exterior.begin(); it != qr_exterior.end(); ++it)
  //     {
  //       const auto& qr = it->second;
  //       const std::size_t cell_index = it->first;
  //       const auto& nn = facet_normals.at(cell_index);

  //       for (std::size_t i = 0; i < qr.second.size(); ++i)
  //       {
  //         const Point p(qr.first[2*i], qr.first[2*i+1]);
  //         std::cout << plot(p,"'k.','markersize',12");
  //         const Point n(nn[2*i], nn[2*i+1]);
  //         const double d = 0.1;
  //         std::cout << drawarrow(p, p+d*n); ///, colors[cutting_cell_part]);
  //       }
  //       std::cout << std::endl;
  //     }
  //   }
      
  // }
  
  //------------------------------------------------------------------------------
  // inline void plot_normals_interface(const MultiMesh& multimesh)
  // {
  //   std::cout << "% "<< __FUNCTION__ << std::endl;
  //   const std::vector<std::string> colors = {{ "'b'", "'g'", "'r'" }};
  //   const std::vector<std::string> marker = {{ "'.'", "'o'", "'x'" }};

  //   for (std::size_t part = 0;
  //        part < std::min<std::size_t>(3, multimesh.num_parts());
  //        part++)
  //     // const std::size_t part = 1;
  //   {
  //     std::cout << "% part " << part << ' ' <<std::endl;
  //     const auto& cmap = multimesh.collision_map_cut_cells(part);
  //     const auto& qr_interface = multimesh.quadrature_rules_interface(part);
  //     const auto& normals = multimesh.facet_normals_interface(part);

  //     for (auto it = cmap.begin(); it != cmap.end(); ++it)
  //     {
  //       const unsigned int cut_cell_index = it->first;
  //       const auto& cutting_cells = it->second;

  //       const Cell cut_cell(*multimesh.part(part), cut_cell_index);
  //       std::cout << drawtriangle(cut_cell, colors[part]);

  //       const auto& qr = multimesh.quadrature_rules_cut_cells(part, cut_cell_index); PPause;
  //       cout_qr(qr, colors[part]);

  //       // Iterate over cutting cells
  //       for (auto jt = cutting_cells.begin(); jt != cutting_cells.end(); jt++)
  //       {
  //         const std::size_t cutting_cell_part = jt->first;

  //         const Cell cutting_cell(*multimesh.part(cutting_cell_part), jt->second);
  //         std::cout << drawtriangle(cutting_cell, colors[cutting_cell_part]);

  //         // Get quadrature rule for interface part defined by
  //         // intersection of the cut and cutting cells
  //         const std::size_t k = jt - cutting_cells.begin();
  //         const auto& qr = qr_interface.at(cut_cell_index)[k];
  //         const auto& nn = normals.at(cut_cell_index)[k];

  //         for (std::size_t i = 0; i < qr.second.size(); ++i)
  //         {
  //           const Point p(qr.first[2*i], qr.first[2*i+1]);
  //           std::cout << plot(p,"'k.','markersize',12");
  //           const Point n(nn[2*i], nn[2*i+1]);
  //           const double d = 0.1;
  //           std::cout << drawarrow(p, p+d*n, colors[cutting_cell_part]);
  //         }
  //       }
  //       std::cout << std::endl;
  //     }

  //     // for (const auto cell_no: multimesh.cut_cells(part))
  //     // {
  //     //   const auto qrmap = multimesh.quadrature_rule_interface(part).find(cell_no);
  //     //   const std::vector<quadrature_rule> qr = qrmap->second;

  //     //   const auto fnmap = multimesh.facet_normals_interface(part).find(cell_no);
  //     //   const std::vector<std::vector<double>> normals = fnmap->second;

  //     //   //std::cout << qr.size() << ' ' << normals.size() << std::endl;
  //     //   dolfin_assert(qr.size() == normals.size());

  //     //   for (std::size_t i = 0; i < qr.size(); ++i)
  //     //   {
  //     // 	for (std::size_t j = 0; j < qr[i].second.size(); ++j)
  //     // 	{
  //     // 	  const Point p(qr[i].first[2*j], qr[i].first[2*j+1]);
  //     // 	  std::cout << plot(p,"'k.'");
  //     // 	  const Point n(normals[i][2*j],normals[i][2*j+1]);
  //     // 	  const double d = 0.01;
  //     // 	  std::cout << drawarrow(p, p+d*n);
  //     // 	}
  //     // 	std::cout << std::endl;
  //     //   }
  //     // }

  //   }
  // }
  //------------------------------------------------------------------------------
  inline void evaluate_at_qr(const MultiMesh& mm,
			     const Expression& uexact,
			     const MultiMeshFunction& uh)
  {
    std::cout << __FUNCTION__ << std::endl;
    double maxee = -1;

    for (std::size_t part = 0; part < mm.num_parts(); ++part)
    {
      std::cout << "\npart " << part << std::endl;

      // get vertex values
      std::vector<double> vertex_values;
      uh.part(part)->compute_vertex_values(vertex_values, *mm.part(part));

      const std::vector<std::string> colors = {{ "'b'", "'g'", "'r'" }};
      std::vector<std::size_t> cells;

      // cells colliding with the cut cells
      const auto collision_map = mm.collision_map_cut_cells(part);

      // loop over cut cells
      for (const auto cut_cell_no: mm.cut_cells(part))
      {
	// all qr on cut_cell_no
	const auto qr = mm.quadrature_rules_cut_cells(part, cut_cell_no); PPause;

	// loop over qr
	for (std::size_t i = 0; i < qr.second.size(); ++i)
	{
	  const Point p(qr.first[2*i], qr.first[2*i+1]);
	  const double uhval = (*uh.part(part))(p);
	  const double uexactval = uexact(p);
	  const double ee = std::abs(uhval-uexactval);
	  maxee = ee > maxee ? ee : maxee;
	  std::cout << p.x()<<' '<<p.y()<<' '<<uhval<<' '<<uexactval<<' '<<ee<<' '<<maxee<<std::endl;

	  //   // if evaluated function big...
	  //   if (std::abs(uhval) > 1)
	  //   {
	  //     // save cell no
	  //     cells.push_back(cut_cell_no);
	  //     const std::string color = qr.second[i] > 0 ? "'.'" : "'x'";
	  //     std::cout << matlabplot(p,color) <<" % " << qr.second[i] << ' '
	  // 	      << /\*std::setprecision(15) <<*\/ uhval << " (";

	  //     // print nodal uh values
	  //     const Cell cell(*mm.part(part), cut_cell_no);
	  //     for (std::size_t j = 0; j < cell.num_vertices(); ++j)
	  //       std::cout << cell.entities(0)[j] << ' '<<vertex_values[cell.entities(0)[j]] <<' ';
	  //     std::cout << ")"<<std::endl;
	  //   }
	}
      }

      // // make cell numbers unique
      // std::sort(cells.begin(), cells.end());
      // const auto new_end = std::unique(cells.begin(), cells.end());
      // cells.erase(new_end, cells.end());

      // // loop over all cells with large uh values
      // for (const auto cell_no: cells)
      // {
      // 	std::cout << "% cell with large uh:"<<std::endl;
      // 	const Cell cell(*mm.part(part), cell_no);
      // 	std::cout << drawtriangle(cell);

      // 	// compute net weight (~visible area)
      // 	const auto qr = mm.quadrature_rule_cut_cell(part, cell_no);
      // 	double net_weight = 0;
      // 	std::cout << " % ";
      // 	for (const auto w: qr.second)
      // 	{
      // 	  net_weight += w;
      // 	  std::cout << ' '<<w;
      // 	}
      // 	std::cout << "\n% net weight = " << net_weight << std::endl;

      // 	// also display all colliding cells
      // 	const auto it = collision_map.find(cell_no);
      // 	dolfin_assert(it->first == cell_no);
      // 	std::cout << "% colliding:"<<std::endl;
      // 	for (const auto cpair: it->second)
      // 	{
      // 	  const Cell cutting_cell(*mm.part(cpair.first), cpair.second);
      // 	  std::cout << drawtriangle(cutting_cell,colors[cpair.first]);
      // 	}
      // }

    }
    std::cout << "max error in qr points " << maxee << std::endl;
    PPause;
  }

  //------------------------------------------------------------------------------
  template<class TFunctionSpace>// eg P1::FunctionSpace
  inline void find_max(const MultiMesh& multimesh,
		       const MultiMeshFunction& u,
		       std::vector<double>& maxvals_parts,
		       std::size_t step = 0
		       // ,
		       // File& uncut0_file, File& uncut1_file, File& uncut2_file,
		       // File& cut0_file, File& cut1_file, File& cut2_file,
		       // File& covered0_file, File& covered1_file, File& covered2_file
		       )

  {
    std::cout << __FUNCTION__ << std::endl;
    std::cout << "\tSolution: max min step " << step <<' ' << u.vector()->max() << ' ' << u.vector()->min() << std::endl;

    maxvals_parts.assign(multimesh.num_parts(), -9e99);

    for (std::size_t part = 0; part < multimesh.num_parts(); ++part)
    {
      // get max on vertex values
      std::vector<double> vertex_values;
      u.part(part)->compute_vertex_values(vertex_values,
					  *multimesh.part(part));
      const double maxvv = *std::max_element(vertex_values.begin(),
					     vertex_values.end());

      // get max on uncut, cut and covered
      const std::vector<std::vector<unsigned int>> cells
	= {{ multimesh.uncut_cells(part),
	     multimesh.cut_cells(part),
	     multimesh.covered_cells(part) }};
      const std::vector<std::string> type = {{ "uncut", "cut", "covered" }};
      std::vector<double> maxvals(cells.size(), 0);

      for (std::size_t k = 0; k < cells.size(); ++k)
      {
	std::cout << "part " << part << " "<<k << ' '<<type[k]<< std::endl;
	if (cells[k].size())
	{
	  // Create meshfunction using markers
	  auto mesh_part = std::make_shared<Mesh>(*multimesh.part(part));
	  auto foo = std::make_shared<MeshFunction<std::size_t> >(mesh_part, mesh_part->topology().dim());
	  foo->set_all(0); // dummy
	  for (const auto cell: cells[k])
	    foo->set_value(cell, k+1);

	  // Create submesh out of meshfunction
	  auto sm = std::make_shared<SubMesh>(*multimesh.part(part), *foo, k+1);

	  // Interpolate on submesh
	  auto V = std::make_shared<TFunctionSpace>(sm);
	  auto usm = std::make_shared<Function>(V);

	  // test
	  usm->set_allow_extrapolation(true);

	  usm->interpolate(*u.part(part));

	  // Get max values on submesh
	  std::vector<double> vertex_values;
	  usm->compute_vertex_values(vertex_values);
	  maxvals[k] = *std::max_element(vertex_values.begin(), vertex_values.end());

	  // if (part == 0)
	  //   if (k == 0 or k == 1) {
	  //     std::cout << k <<std::endl;
	  //     for (const auto cell: cells[k])
	  // 	std::cout << cell << ' ';
	  //     std::cout << std::endl;
	  //   }

	  // if (marker == 1 and part == 0) {
	  //   for (const auto v: vertex_values)
	  //     std::cout << v<<' ';
	  //   std::cout << std::endl;
	  // }

	  // // save
	  // switch(k) {
	  // case 0: { // uncut
	  //   if (part == 0) uncut0_file << (*usm);
	  //   else if (part == 1) uncut1_file << (*usm);
	  //   else if (part == 2) uncut2_file << (*usm);
	  //   break;
	  // }
	  // case 1: { // cut
	  //   if (part == 0) cut0_file << (*usm);
	  //   else if (part == 1) cut1_file << (*usm);
	  //   else if (part == 2) cut2_file << (*usm);
	  //   break;
	  // }
	  // case 2: { // covered
	  //   if (part == 0) covered0_file << (*usm);
	  //   else if (part == 1) covered1_file << (*usm);
	  //   else if (part == 2) covered2_file << (*usm);
	  // }
	  // }
	}
      }

      std::cout << "\tpart " << part
		<< " step " << step
		<< " all vertices " << maxvv
		<< " uncut " << maxvals[0]
		<< " cut " << maxvals[1]
		<< " covered " << maxvals[2] << std::endl;

      maxvals_parts[part] = std::max(std::max(maxvals[0], maxvals[1]), maxvals[2]);
      //if (maxvals[0] < 1) { exit(0); }
    }

  }

  //------------------------------------------------------------------------------
  template<class T>
  inline std::string write_matrix_raw(const T& A)
  {
    std::stringstream s;
    for (std::size_t r = 0; r < A.size(0); ++r)
    {
      std::vector<std::size_t> columns;
      std::vector<double> values;
      A.getrow(r, columns, values);
      for (std::size_t c = 0; c < columns.size(); ++c)
	s << r << ' '<< columns[c] << ' ' << values[c] << std::endl;
    }
    return s.str();
  }

  //-----------------------------------------------------------------------------
  template<class T>
  inline void write_matrix(const std::string& filename,
			   const T& A)
  {
    std::ofstream f(filename);
    if (!f.good()) { std::cout << "trying to open file " << filename << " from function " << __PRETTY_FUNCTION__ << " but it failed\n"; exit(1); }
    f << tools::write_matrix_raw(A);
    f.close();
  }

  //-----------------------------------------------------------------------------
  inline void check_finite(const Matrix& A)

  {
    std::vector<std::size_t> cols;
    std::vector<double> vals;
    for (std::size_t r = 0; r < A.size(0); ++r)
    {
      A.getrow(r, cols, vals);
      for (std::size_t i = 0; i < cols.size(); ++i)
	if (!std::isfinite(vals[i]))
	  std::cout << r << ' ' << cols[i] << ' ' << vals[i] << std::endl;
    }
  }

  //-----------------------------------------------------------------------------
  inline std::string generate_test(dolfin::Point p0,
				   dolfin::Point p1,
				   dolfin::Point q0,
				   dolfin::Point q1,
				   std::string function)
  {
    std::stringstream ss;
    ss.precision(std::numeric_limits<long double>::digits10+1);
    ss << "@skip_in_parallel\n"
       << "def test_segment_segment_XX():\n"
       << "    \"Case that fails CGAL comparison. We get a different intersection point but still correct area.\"\n"
       << "    p0 = Point(" << p0[0]<<','<<p0[1]<<")\n"
       << "    p1 = Point(" << p1[0]<<','<<p1[1]<<")\n"
       << "    q0 = Point(" << q0[0]<<','<<q0[1]<<")\n"
       << "    q1 = Point(" << q1[0]<<','<<q1[1]<<")\n"
       << "    intersection = IntersectionConstruction.intersection_segment_segment_2d(p0, p1, q0, q1)\n\n"
       << "    # The intersection should according to CGAL be\n"
       << "    cgal = Point(\n\n"
       << "    # We get\n"
       << "    computed = Point(\n\n"
       << "    assert len(intersection) == 1\n"
       << "    assert (abs(intersection[0][0] - cgal[0]) < DOLFIN_EPS and abs(intersection[0][1] - cgal[1]) < DOLFIN_EPS) or \n"
       << "        (abs(intersection[0][0] - computed[0]) < DOLFIN_EPS and abs(intersection[0][1] - computed[1]) < DOLFIN_EPS)\n";

    return ss.str();
  }
  /* //----------------------------------------------------------------------------- */
  /* inline double condition_number(const std::size_t step, */
  /* 				 const std::shared_ptr<PETScMatrix> A) */
  /* { */
  /*   std::cout << "Condition number using PETSc/SLEPc" << std::endl; */
  /*   // See example http://slepc.upv.es/documentation/current/src/svd/examples/tutorials/ex8.c.html */

  /*   SVD svd; */
  /*   PetscInt nconv1, nconv2; */
  /*   PetscReal sigma_1, sigma_n; */

  /*   // Setup */
  /*   PETScOptions::set("svd_view"); */
  /*   PETScOptions::set("svd_monitor"); */
  /*   //PETScOptions::set<int>("svd_ncv", 500); */
  /*   PETScOptions().set("svd_type", "trlanczos"); */
  /*   // PETScOptions().set("svd_eps_type", "krylovschur"); */
  /*   // PETScOptions().set("svd_monitor_draw", ""); */

  /*   SVDCreate(PETSC_COMM_WORLD,&svd); */
  /*   SVDSetOperator(svd, A->mat()); */
  /*   SVDSetFromOptions(svd); */

  /*   nconv2 = -1; */
  /*   PetscInt svd_ncv = 16; */
  /*   PetscReal tol = 1e-3; */
  /*   PetscInt nit = 1000; */
  /*   SVDSetTolerances(svd,tol,nit); */

  /*   while (nconv2 <= 0) */
  /*   { */
  /*     std::cout << nconv2 << ' '<< svd_ncv << std::endl; */
  /*     SVDSetDimensions(svd,1,svd_ncv,PETSC_DEFAULT); */

  /*     // Solve the svd problem for the smallest singular value */
  /*     SVDSetWhichSingularTriplets(svd,SVD_SMALLEST); */
  /*     SVDSolve(svd); */
  /*     SVDGetConverged(svd,&nconv2); */
  /*     if (nconv2 > 0) { */
  /* 	SVDGetSingularTriplet(svd,0,&sigma_n,NULL,NULL); */
  /*     } else { */
  /* 	PetscPrintf(PETSC_COMM_WORLD," Unable to compute small singular value!\n\n"); */
  /* 	svd_ncv *= 2; */
  /*     } */
  /*   } */

  /*   // Find largest */
  /*   SVDSetWhichSingularTriplets(svd,SVD_LARGEST); */
  /*   SVDSolve(svd); */
  /*   SVDGetConverged(svd,&nconv1); */
  /*   if (nconv1 > 0) { */
  /*     SVDGetSingularTriplet(svd,0,&sigma_1,NULL,NULL); */
  /*   } else { */
  /*     PetscPrintf(PETSC_COMM_WORLD," Unable to compute large singular value!\n\n"); */
  /*     exit(1); */
  /*   } */

  /*   SVDDestroy(&svd); */

  /*   const double condno = (double)(sigma_1/sigma_n); */
  /*   std::cout << "singular values " << sigma_n<<' '<<sigma_1 << "  "<<condno << std::endl; */

  /*   return condno; */
  /* } */

  //-----------------------------------------------------------------------------
  inline std::string print_std_vector(const std::vector<Point>& p,
				      const std::string& var)
  {
    std::stringstream ss;
    //ss.precision(16);
    ss << std::setprecision(std::numeric_limits<long double>::digits10+2);
    ss << "const std::vector<Point> " << var << " = {{";
    for (std::size_t i = 0; i < p.size()-1; ++i)
    {
      ss << "Point(" << p[i].x() <<',' << p[i].y() << ',' << p[i].z() << "),\n";
    }
    std::size_t i = p.size()-1;
    ss << "Point(" << p[i].x() <<',' << p[i].y() << ',' << p[i].z() << ") }};\n";
    return ss.str();
  }
  //-----------------------------------------------------------------------------
  inline std::string print_std_vector(const Cell& cell,
				      const std::string& var)
  {
    return print_std_vector(convert(cell), var);
  }
  //-----------------------------------------------------------------------------

}

#endif
