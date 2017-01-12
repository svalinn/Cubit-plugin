#include "mcnp2cad.hpp"
#include "DAGMCExportCommand.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <cctype>
#include <vector>
#include <set>
#include <map>
#include <sstream>
#include <algorithm>

#include <cassert>

//MCNP2CAD includes
#include "mcnp2cad/geometry.hpp"
#include "mcnp2cad/MCNPInput.hpp"
#include "mcnp2cad/options.hpp"
#include "mcnp2cad/volumes.hpp"
#include "mcnp2cad/ProgOptions.hpp"
#include "mcnp2cad/version.hpp"
#include "CubitInterface.hpp"

// CGM includes
#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "ModelQueryEngine.hpp"
#include "GMem.hpp"

#include "RefEntityName.hpp"
#include "RefEntity.hpp"

#include "Body.hpp"
#include "Surface.hpp"
#include "Curve.hpp"

#include "RefGroup.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"

#include "SenseEntity.hpp"

// MOAB includes[
#include "MBTagConventions.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/GeomTopoTool.hpp"


typedef std::vector<iBase_EntityHandle> entity_collection_t;

static bool boundBoxesIntersect( iBase_EntityHandle h1, iBase_EntityHandle h2 ){

  Vector3d h1_min, h1_max, h2_min, h2_max;
  int igm_result;
  
  iGeom_getEntBoundBox( h1, h1_min.v, h1_min.v+1, h1_min.v+2, h1_max.v, h1_max.v+1, h1_max.v+2, &igm_result );
  CHECK_IGEOM( igm_result, "Getting bounding box h1" );
  iGeom_getEntBoundBox( h2, h2_min.v, h2_min.v+1, h2_min.v+2, h2_max.v, h2_max.v+1, h2_max.v+2, &igm_result );
  CHECK_IGEOM( igm_result, "Getting bounding box h2" );

  bool ret = false;

  for( int i = 0; i < 3 && ret == false; ++i ){
    ret = ret || ( h1_min.v[i] > h2_max.v[i] );
    ret = ret || ( h2_min.v[i] > h1_max.v[i] );
  }

  return !ret;

}

static bool intersectIfPossible( iBase_EntityHandle h1, iBase_EntityHandle h2, iBase_EntityHandle* result, 
                                 bool delete_on_failure = true)
{
//TODO Find some way to determine successful intersection.
  int igm_result;

  iGeom_intersectEnts( h1, h2, result, &igm_result);
  
  if( igm_result == iBase_SUCCESS ){
    return true;
  }
  else{
    if( delete_on_failure ){
      iGeom_deleteEnt( h1, &igm_result);
      CHECK_IGEOM(igm_result, "deleting an intersection candidate");
      iGeom_deleteEnt( h2, &igm_result);
      CHECK_IGEOM(igm_result, "deleting an intersection candidate");
    }
    return false;
  }
}

class MCNP2CAD::GeometryContext {

  /** 
   * Metadata and naming: 
   * The NamedGroup and NamedEntity mappings are used to keep track of metadata
   * on particular entity handles that map to MCNP cells.  
   * EntityHandles change frequently as CSG operations are performed on volumes,
   * so these mappings must be updated, by calling updateMaps(), whenever an 
   * EntityHandle changes.
   */

protected:

  // note: this appears slow, since it's called for all cells and constructs a string
  // for each.  A lookup table would probably be faster, but there are never more than
  // a few thousand cells.
  std::string materialName( int mat, double rho ){
    std::string ret;
    std::stringstream formatter;
//    if(Gopt.uwuw_names){
      bool mass_density = false;
      if (rho <= 0){
        mass_density = true;
        rho = -rho;
      }
      char rho_formatted [50];
      sprintf(rho_formatted, "%E", rho);
      formatter << "mat:m" << mat;
      if(mass_density)
          formatter << "/rho:" <<rho_formatted;
      else
          formatter << "/atom:" << rho_formatted;
//    }
//    else
//      formatter << "mat_" << mat << "_rho_" << rho;
    formatter >> ret;
    return ret;
  }

  std::string importanceName( char impchar, double imp ){
    std::string ret; 
    std::stringstream formatter;
    formatter << "imp." << impchar << "_" << imp;
    formatter >> ret;
    return ret;
  }

  class NamedGroup {
  protected:
    std::string name;
    entity_collection_t entities;
  public:
    NamedGroup( ) : name("") {}
    NamedGroup( std::string name_p ):
      name(name_p)
    {}

    const std::string& getName() const { return name; }
    const entity_collection_t& getEntities() const { return entities; }

    void add( iBase_EntityHandle new_handle ){
      entities.push_back(new_handle);
    }

    void update( iBase_EntityHandle old_h, iBase_EntityHandle new_h ){
      entity_collection_t::iterator i = std::find( entities.begin(), entities.end(), old_h );
      if( i != entities.end() ){
        if( new_h ){
          *i = new_h;
        }
        else {
          entities.erase( i );
        }
      }
    }

    bool contains( iBase_EntityHandle handle ) const {
      return std::find( entities.begin(), entities.end(), handle ) != entities.end(); 
    }

  };

  class NamedEntity { 
  protected:
    iBase_EntityHandle handle;
    std::string name;
  public:
    NamedEntity( iBase_EntityHandle handle_p, std::string name_p = "" ):
      handle(handle_p), name(name_p)
    {}
    virtual ~NamedEntity(){}

    const std::string& getName() const { return name; }
    iBase_EntityHandle getHandle() const{ return handle; }

    void setHandle( iBase_EntityHandle new_h ) {
      handle = new_h;
    }

    static NamedEntity* makeCellIDName( iBase_EntityHandle h, int ident ){
      NamedEntity* e = new NamedEntity(h);
      std::stringstream formatter;
      formatter << "MCNP_ID_" << ident;
      formatter >> e->name;
      return e;
    }

  };

protected:
  iGeom_Instance& igm;
  InputDeck& deck;
  double world_size;
  int universe_depth;

  std::map< std::string, NamedGroup* > named_groups;
  std::vector< NamedEntity* > named_cells;

  NamedGroup* getNamedGroup( const std::string& name ){
    if( named_groups.find( name ) == named_groups.end() ){
      named_groups[ name ] = new NamedGroup( name );
//      if( OPT_DEBUG ) std::cout << "New named group: " << name 
 //                               << "num groups now " << named_groups.size() << std::endl;
    }
    return named_groups[ name ];
  }

public:
  GeometryContext( iGeom_Instance& igm_p, InputDeck& deck_p ) :
    igm(igm_p), deck(deck_p), world_size(0.0), universe_depth(0)
  {}

  bool defineLatticeNode( CellCard& cell, iBase_EntityHandle cell_shell, iBase_EntityHandle lattice_shell,
                          int x, int y, int z, entity_collection_t& accum );
  

  entity_collection_t defineCell( CellCard& cell, bool defineEmbedded, iBase_EntityHandle lattice_shell );
  entity_collection_t populateCell( CellCard& cell, iBase_EntityHandle cell_shell, iBase_EntityHandle lattice_shell );
 
  //TODO Change this; don't have initial values in declaration.
  entity_collection_t defineUniverse( int universe, iBase_EntityHandle container, const Transform* transform = NULL );
  

  void addToVolumeGroup( iBase_EntityHandle cell, const std::string& groupname );
  void setVolumeCellID( iBase_EntityHandle cell, int ident);
  void setMaterial( iBase_EntityHandle cell, int material, double density ){
 //   if( Gopt.tag_materials ){
        addToVolumeGroup( cell, materialName(material,density) );
  //  }
  }

  void setImportances( iBase_EntityHandle cell, const std::map<char, double>& imps ){
//    if( Gopt.tag_importances ){
      for( std::map<char, double>::const_iterator i = imps.begin();
           i != imps.end(); ++i )
      {
          char impchar = (*i).first;
          double imp = (*i).second;
          addToVolumeGroup( cell, importanceName( impchar, imp ) );
      }
//    }
  }

  void createGeometry( ){


    int igm_result;

    InputDeck::cell_card_list    cells     = deck.getCells();
    InputDeck::surface_card_list surfaces  = deck.getSurfaces();
    InputDeck::data_card_list    datacards = deck.getDataCards();
    int k = 0;
    int ident, facet;

    // estimate how large the geometry will need to be to accomodate all the surfaces
    for( InputDeck::surface_card_list::iterator i = surfaces.begin(); i!=surfaces.end(); ++i){
      // catch all exceptions from makeSurface at this time; if they exist, they will
      // more properly be displayed to the user at a later time.  Right now we just want
      // to estimate a size and failures can be ignored.
      try{
        if( surfaces.at(k)->getIdent() < 0 ){
          ident = surfaces.at(k)->getIdent()/10;
          facet = -( surfaces.at(k)->getIdent() - 10*ident );
        }
        else{
          facet = 0;
        }
        k++;

        world_size = std::max( world_size, makeSurface( *i, NULL, facet ).getFarthestExtentFromOrigin() );
      } catch(std::runtime_error& e){}
    }

    // translations can extend the size of the world
    double translation_addition = 0;
    for( InputDeck::data_card_list::iterator i = datacards.begin(); i!=datacards.end(); ++i){
      DataCard* c = *i;
      if( c->getKind() == DataCard::TR ){
        double tform_len = dynamic_cast<DataRef<Transform>*>(c)->getData().getTranslation().length();
        translation_addition = std::max (translation_addition, tform_len );
      }
    }

    for( InputDeck::cell_card_list::iterator i = cells.begin(); i!=cells.end(); ++i){
      CellCard* c = *i;
      // translations can come from TRCL data
      if( c->getTrcl().hasData() ){
        double tform_len = c->getTrcl().getData().getTranslation().length();
        translation_addition = std::max( translation_addition, tform_len );
      }
      // translations can also come from fill nodes.  This implementation does *not* take
      // lattices into account, as they are assumed to be contained within other volumes.
      if( c->hasFill() && c->getFill().getOriginNode().hasTransform() ){
        double tform_len = c->getFill().getOriginNode().getTransform().getTranslation().length();
        translation_addition = std::max( translation_addition, tform_len );
      }
    }

    world_size += translation_addition;
    world_size *= 1.2;

    std::cout << "World size: " << world_size << " (trs added " << translation_addition << ")" << std::endl;

    iBase_EntityHandle graveyard = NULL, graveyard_boundary = NULL;
//    if( Gopt.make_graveyard ){
      graveyard = createGraveyard ( graveyard_boundary ); 
//    }

    std::cout << "Defining geometry..." << std::endl;

    entity_collection_t defined_cells = defineUniverse( 0, graveyard_boundary );
    if( graveyard ){ defined_cells.push_back(graveyard); }

    size_t count = defined_cells.size();
    iBase_EntityHandle *cell_array = new iBase_EntityHandle[ count ];
    for( unsigned int i = 0; i < count; ++i ){
      cell_array[i] = defined_cells[i];
    }

    // if( OPT_DEBUG ){ mapSanityCheck(cell_array, count); }
    tagGroups();

    //  if( Gopt.tag_cell_IDs ){
    tagCellIDsAsEntNames();
    //  }


    //  if ( Gopt.imprint_geom ) {
    std::cout << "Imprinting all...\t\t\t" << std::flush;
    iGeom_imprintEnts( cell_array, count, &igm_result );
    CHECK_IGEOM( igm_result, "Imprinting all cells" );
    std::cout << " done." << std::endl;

    double tolerance = world_size / 1.0e7;
    //    if( Gopt.override_tolerance ){
    //      tolerance = Gopt.specific_tolerance;
    //    }

    //    if ( Gopt.merge_geom ) {
    std::cout << "Merging, tolerance=" << tolerance << "...\t\t" << std::flush;
    iGeom_mergeEnts( cell_array, count,  tolerance, &igm_result );
    CHECK_IGEOM( igm_result, "Merging all cells" );
    std::cout << " done." << std::endl;
    //    }
    //  }


    /*
       std::string outName = Gopt.output_file;
       std::cout << "Saving file \"" << outName << "\"...\t\t\t" << std::flush;
       iGeom_save( outName.c_str(), ""/*, &igm_result* /, outName.length(), 0 );
       CHECK_IGEOM( igm_result, "saving the output file "+outName );
       std::cout << " done." << std::endl;
     */

  }

  void updateMaps ( iBase_EntityHandle old_cell, iBase_EntityHandle new_cell );

  bool mapSanityCheck( iBase_EntityHandle* cells, size_t count );

  void tagGroups( );
  void tagCellIDsAsEntNames();

  std::string uprefix() { 
    return std::string( universe_depth, ' ' );
  }

  iBase_EntityHandle createGraveyard( iBase_EntityHandle& boundary );
//  void createGeometry( );

};

entity_collection_t MCNP2CAD::GeometryContext::defineUniverse( int universe, iBase_EntityHandle container = NULL, 
                                                     const Transform* transform /*= NULL*/ )
{

//  if( OPT_VERBOSE ) std::cout << uprefix() << "Defining universe " << universe << std::endl;
  universe_depth++;

  InputDeck::cell_card_list u_cells = deck.getCellsOfUniverse( universe );
  entity_collection_t subcells;

  iBase_EntityHandle lattice_shell = NULL;
  if( u_cells.size() == 1 && u_cells[0]->isLattice() ){
    lattice_shell = container;
    // reverse-transform the containing volume before using it as a lattice boundary
    if(transform){
      lattice_shell = applyReverseTransform( *transform, lattice_shell );
    }
  }

  // define all the cells of this universe
  for( InputDeck::cell_card_list::iterator i = u_cells.begin(); i!=u_cells.end(); ++i){
    entity_collection_t tmp = defineCell( *(*i), true, lattice_shell );
    for( size_t i = 0; i < tmp.size(); ++i){
      subcells.push_back( tmp[i] );
    }
  }
  
  if( transform ){
    for( size_t i = 0; i < subcells.size(); ++i){
      subcells[i] = applyTransform( *transform, subcells[i] );      
    }
  }

  if( container && !lattice_shell ){
    
    int igm_result;
    
    for( size_t i = 0; i < subcells.size(); ++i ){
      
//      if( OPT_DEBUG ) std::cout << uprefix() << "Bounding a universe cell..." << std::flush;    
      bool subcell_removed = false;

      if( boundBoxesIntersect( subcells[i], container )){
        iBase_EntityHandle container_copy;
        iGeom_copyEnt( container, &container_copy, &igm_result);
        CHECK_IGEOM( igm_result, "Copying a universe-bounding cell" );
        
        iBase_EntityHandle subcell_bounded;
        bool valid_result = intersectIfPossible( container_copy, subcells[i], &subcell_bounded );
        if( valid_result ){
          updateMaps( subcells[i], subcell_bounded );
          subcells[i] = subcell_bounded;
//          if( OPT_DEBUG ) std::cout << " ok." <<  std::endl;
        }
        else{
          subcell_removed = true;
        }

      }
      else{
        // bounding boxes didn't intersect, delete subcells[i].
        // this suggests invalid geometry, but we can continue anyway.
        iGeom_deleteEnt( subcells[i], &igm_result );
        CHECK_IGEOM( igm_result, "Deleting a subcell that didn't intersect a parent's bounding box (strange!)" );
        subcell_removed = true;
      }
      
      if( subcell_removed ){
        updateMaps( subcells[i], NULL );
        subcells.erase( subcells.begin()+i );
        i--;
//        if( OPT_DEBUG ) std::cout << " removed." << std::endl;
      }
      
    }
        
    iGeom_deleteEnt( container, &igm_result );
    CHECK_IGEOM( igm_result, "Deleting a bounding cell" );
  }

  universe_depth--;
//  if( OPT_VERBOSE ) std::cout << uprefix() << "Done defining universe " << universe << std::endl;

  return subcells;
 
}

iBase_EntityHandle MCNP2CAD::GeometryContext::createGraveyard( iBase_EntityHandle& inner_copy ) {
  iBase_EntityHandle inner, outer, graveyard;
  int igm_result;

  double inner_size = 2.0 * world_size;
  iGeom_createBrick( inner_size, inner_size, inner_size, &inner, &igm_result );
  CHECK_IGEOM( igm_result, "Making graveyard" );
  
  iGeom_copyEnt( inner, &inner_copy, &igm_result );
  CHECK_IGEOM( igm_result, "Copying graveyard" );
  
  double outer_size = 2.0 * ( world_size + (world_size / 50.0) );
  iGeom_createBrick( outer_size, outer_size, outer_size, &outer, &igm_result );
  CHECK_IGEOM( igm_result, "Making outer graveyard" );

  iGeom_subtractEnts( outer, inner, &graveyard, &igm_result );
  CHECK_IGEOM( igm_result, "subtracting graveyard" );
  
  addToVolumeGroup( graveyard, "graveyard" );
  
  // reset world size to a sphere that bounds the inner shell of this graveyard
  world_size *= sqrt(3.0);
//if( OPT_DEBUG ) std::cout << "Spherical world size for graveyard: " << world_size << std::endl;
  
  return graveyard;

}


void MCNP2CAD::GeometryContext::addToVolumeGroup( iBase_EntityHandle cell, const std::string& name ){

  NamedGroup* group = getNamedGroup( name );
  group->add( cell );

//  if( OPT_DEBUG ){ std::cout << uprefix() 
//                             << "Added cell to volgroup " << group->getName() << std::endl; }
}

entity_collection_t MCNP2CAD::GeometryContext::defineCell(  CellCard& cell,  bool defineEmbedded = true, 
                                                  iBase_EntityHandle lattice_shell = NULL )
{
  int ident = cell.getIdent();
  const CellCard::geom_list_t& geom = cell.getGeom();
 
//  if( OPT_VERBOSE ) std::cout << uprefix() << "Defining cell " << ident << std::endl;

  int igm_result;

  entity_collection_t tmp;

  std::vector<iBase_EntityHandle> stack;
  for(CellCard::geom_list_t::const_iterator i = geom.begin(); i!=geom.end(); ++i){
    
    const CellCard::geom_list_entry_t& token = (*i);
    switch(token.first){
    case CellCard::CELLNUM:
      // a cell number appears in a geometry list only because it is being complemented with the # operator
      // thus, when defineCell is called on it, set defineEmbedded to false
      tmp = defineCell( *(deck.lookup_cell_card(token.second)), false);
      assert(tmp.size() == 1);
      stack.push_back( tmp.at(0) );
      break;
    case CellCard::SURFNUM:
      {      
        int surface = token.second;
        bool positive = true;
        if( surface < 0){
          positive = false; surface = -surface;
        }
        try{
          SurfaceVolume& surf = makeSurface( deck.lookup_surface_card( surface ) );
          iBase_EntityHandle surf_handle = surf.define( positive, world_size );
          stack.push_back(surf_handle);
        }
        catch(std::runtime_error& e) { std::cerr << e.what() << std::endl; }
      }
      break;
    case CellCard::MBODYFACET:
      {
        int identifier = -std::abs( token.second );
        int surfacenum = -identifier / 10;
        int facet = -identifier - ( surfacenum * 10 );

        try{
          SurfaceVolume& surf = makeSurface( deck.lookup_surface_card( identifier ) );
          const std::string& mnemonic = deck.lookup_surface_card( identifier )->getMnemonic();
          bool positive = true;
          if( mnemonic == "rcc" || mnemonic == "rec" ){
            if( ( token.second < 0 ) ^ ( facet == 3 ) ){
              positive = false;
            }
          }
          else if( mnemonic == "box" || mnemonic == "rpp" ){
            if( ( token.second < 0 ) ^ ( facet == 2 || facet == 4 || facet == 6 ) ){
              positive = false;
            }
          }
          else if( mnemonic == "hex" || mnemonic == "rhp" ){
            if( ( token.second < 0 ) ^ ( facet == 8 ) ){
              positive = false;
            }
          }

          iBase_EntityHandle surf_handle = surf.define ( positive, world_size );
          stack.push_back(surf_handle);
        }
        catch(std::runtime_error& e) { std::cerr << e.what() << std::endl; }


      }
      break;
    case CellCard::INTERSECT:
      {
        assert( stack.size() >= 2 );
        iBase_EntityHandle s1 = stack.back(); stack.pop_back();
        iBase_EntityHandle s2 = stack.back(); stack.pop_back();
        iBase_EntityHandle result;
        if( intersectIfPossible( s1, s2, &result ) ){
          stack.push_back(result);
        }
        else{
          std::cout << "FAILED INTERSECTION CELL #" << cell.getIdent() << std::endl;
          throw std::runtime_error("Intersection failed");
        }
      }
      break;
    case CellCard::UNION:
      { 
        assert( stack.size() >= 2 );
        iBase_EntityHandle s[2];
        s[0] = stack.back(); stack.pop_back();
        s[1] = stack.back(); stack.pop_back();
        iBase_EntityHandle result;
        iGeom_uniteEnts( s, 2, &result, &igm_result);
        CHECK_IGEOM( igm_result, "Uniting two entities" );
        stack.push_back(result);
      }
      break;
    case CellCard::COMPLEMENT:
      {
        assert (stack.size() >= 1 );
        iBase_EntityHandle world_sphere = makeWorldSphere( world_size);
        iBase_EntityHandle s = stack.back(); stack.pop_back();
        iBase_EntityHandle result;

        iGeom_subtractEnts(  world_sphere, s, &result, &igm_result);
        CHECK_IGEOM( igm_result, "Complementing an entity" );
        stack.push_back(result);
      }
      break;
    default:
      throw std::runtime_error( "Unexpected token while evaluating cell geometry");
      break;
    }
  }

  assert( stack.size() == 1);

  iBase_EntityHandle cellHandle = stack[0];

  if( cell.getTrcl().hasData() ){
    cellHandle = applyTransform( cell.getTrcl().getData(), cellHandle );
  }

  if( defineEmbedded ){
    return populateCell( cell, cellHandle, lattice_shell );
  }
  else{
    return entity_collection_t( 1, cellHandle );
  }
  
}

typedef struct{ int v[3]; } int_triple;

static std::vector<int_triple> makeGridShellOfRadius( int r, int dimensions ){
  if( r == 0 ){ 
    int_triple v; v.v[0] = v.v[1] = v.v[2] = 0;
    return std::vector<int_triple>(1,v);
  }
  else{
    std::vector<int_triple> ret;
    
    int jmin = dimensions > 1 ? -r : 0;
    int jmax = dimensions > 1 ?  r : 0;
    int kmin = dimensions > 2 ? -r : 0;
    int kmax = dimensions > 2 ?  r : 0;
    for( int i = -r; i <= r; ++i ){
      for( int j = jmin;j <= jmax; ++j ){
        for( int k = kmin; k <= kmax; ++k ){
          if( i == -r || i == r ||
              j == -r || j == r ||
              k == -r || k == r ){
            int_triple v;
            v.v[0] = i; 
            v.v[1] = j;
            v.v[2] = k;
            ret.push_back(v);
          }
        }
      }
    }
    return ret;
  }
}

entity_collection_t MCNP2CAD::GeometryContext::populateCell( CellCard& cell,  iBase_EntityHandle cell_shell, 
                                                   iBase_EntityHandle lattice_shell = NULL )
{
  
//  if( OPT_DEBUG ) std::cout << uprefix() << "Populating cell " << cell.getIdent() << std::endl;


  if( !cell.hasFill() && !cell.isLattice() ){
    // no further geometry inside this cell; set its material
    setVolumeCellID(cell_shell, cell.getIdent());
    if( cell.getMat() != 0 ){ setMaterial( cell_shell, cell.getMat(), cell.getRho() ); }
    if( cell.getImportances().size() ){ setImportances( cell_shell, cell.getImportances()); }
    return entity_collection_t(1, cell_shell );
  }
  else if(cell.hasFill() && !cell.isLattice()){
    // define a simple (non-lattice) fill
    
    const FillNode& n = cell.getFill().getOriginNode();
    int filling_universe = n.getFillingUniverse();
    
/*    if( OPT_DEBUG ){
      std::cout << uprefix() << "Creating cell " << cell.getIdent() 
                << ", which is filled with universe " << filling_universe << std::endl;
    }
    */
    
    // the contained universe is transformed by the FillNode's transform, if any, or
    // else by the cell's TRCL value, if any.
    const Transform* t;
    if( n.hasTransform() ){
      t = &(n.getTransform());
    } else if( cell.getTrcl().hasData() ){
      t = &(cell.getTrcl().getData() );
    } else { 
      t = NULL; 
    }

//    if( OPT_DEBUG && t ) std::cout << uprefix() << " ... and has transform: " << *t << std::endl;

    entity_collection_t subcells = defineUniverse(  filling_universe, cell_shell, t );
 
    return subcells;
     
  }
  else {
    // cell is a lattice, bounded by lattice_shell.  cell_shell is the origin element of the lattice and
    // cell->getLattice() has the lattice parameters.

    assert(lattice_shell);
    
//    if( OPT_VERBOSE ) std::cout << uprefix() << "Creating cell " << cell.getIdent() << "'s lattice" << std::endl;

    entity_collection_t subcells;
        
    const Lattice& lattice = cell.getLattice();
    int num_dims = lattice.numFiniteDirections();
    
//    if( OPT_DEBUG ) std::cout << uprefix() << "  lattice num dims: " << num_dims << std::endl;

    if( lattice.isFixedSize() ){

//      if( OPT_DEBUG ) std::cout << uprefix() << "Defining fixed lattice" << std::endl;

      irange xrange = lattice.getXRange(), yrange = lattice.getYRange(), zrange = lattice.getZRange();

      for( int i = xrange.first; i <= xrange.second; ++i){
        for( int j = yrange.first; j <= yrange.second; ++j ){
          for( int k = zrange.first; k <= zrange.second; ++k ){

//            if( OPT_DEBUG ) std::cout << uprefix() << "Defining lattice node " << i << ", " << j << ", " << k << std::endl;

            /* bool success = */ defineLatticeNode( cell, cell_shell, lattice_shell, i, j, k, subcells );

            if( num_dims < 3 ) break; // from z loop
          }
          if( num_dims < 2 ) break; // from y loop
        }
      }

    }
    else{

//      if( OPT_DEBUG ) std::cout << uprefix() << "Defining infinite lattice" << std::endl;
//      if( OPT_VERBOSE && Gopt.infinite_lattice_extra_effort ) 
//        std::cout << uprefix() << "Infinite lattice extra effort enabled." << std::endl;

      // when extra effort is enabled, initialize done_one to false;
      // the code will keep trying to create lattice elements until at least one 
      // element has been successfully created.
      bool done = false, done_one = false;//!Gopt.infinite_lattice_extra_effort;
      int radius = 0;

      while( !done ){
        
        done = done_one;
        std::vector<int_triple> shell = makeGridShellOfRadius(radius++, num_dims);

        for( std::vector<int_triple>::iterator i = shell.begin(); i!=shell.end(); ++i){
          int x = (*i).v[0];
          int y = (*i).v[1];
          int z = (*i).v[2];
          
//          if( OPT_DEBUG ) std::cout << uprefix() << "Defining lattice node " << x << ", " << y << ", " << z << std::endl;

          bool success = defineLatticeNode( cell, cell_shell, lattice_shell, x, y, z, subcells );
          if( success ){
            done = false;
            done_one = true;
          }

        }       
      }
    }

    int igm_result;
    iGeom_deleteEnt( cell_shell, &igm_result );
    CHECK_IGEOM( igm_result, "Deleting cell shell after building lattice" );
    iGeom_deleteEnt( lattice_shell, &igm_result );
    CHECK_IGEOM( igm_result, "Deleting lattice shell after building lattice" );

    return subcells;
  }
}

void MCNP2CAD::GeometryContext::setVolumeCellID( iBase_EntityHandle cell, int ident ){

  named_cells.push_back( NamedEntity::makeCellIDName(cell, ident) );

}

void MCNP2CAD::GeometryContext::updateMaps( iBase_EntityHandle old_cell, iBase_EntityHandle new_cell ){

  /* update named_groups.  handling of new_cell == NULL case is performed within NamedGroup class */
  for( std::map<std::string,NamedGroup*>::iterator i = named_groups.begin();
       i != named_groups.end(); ++i )
  {
    NamedGroup* group = (*i).second;
    group->update( old_cell, new_cell );
  }

  /* update named entities.*/
  if( new_cell != NULL ){
    for( std::vector< NamedEntity* >::iterator i = named_cells.begin(); 
         i != named_cells.end(); ++i )
      {
        NamedEntity* ne = *i;
        if( ne->getHandle() == old_cell ){
          ne->setHandle( new_cell );
        }
      }
  }
  else{ /* new_cell == NULL (i.e. cell has disappeared) */

    // this case is expected to be uncommon in most geometries, so the fact that erasing
    // from a std::vector is slow should not be a problem.
    std::vector< NamedEntity* >::iterator i = named_cells.begin();
    while( i != named_cells.end() ){
      if( (*i)->getHandle() == old_cell ){
        delete (*i);
        named_cells.erase(i);
      }
      else{
        ++i;
      }
    }
  }

}

void MCNP2CAD::GeometryContext::tagGroups( ){

  // the NamedGroup system used to be used solely to create groups for material specification,
  // but it has since been expanded for importance groups.  Some of the old messages that
  // talk about material groups could be confusing.
  int igm_result;
  
  std::string name_tag_id = "NAME";
  int name_tag_maxlength = 64;
  iBase_TagHandle name_tag;

  iGeom_getTagHandle( igm, name_tag_id.c_str(), &name_tag, &igm_result, name_tag_id.length() );
  CHECK_IGEOM( igm_result, "Looking up NAME tag" );
  
  iGeom_getTagSizeBytes( igm, name_tag, &name_tag_maxlength, &igm_result );
  CHECK_IGEOM( igm_result, "Querying NAME tag length" );
//  if( OPT_DEBUG ) std::cout << "Name tag length: " << name_tag_maxlength << " actual id " << name_tag << std::endl;

  for( std::map<std::string,NamedGroup*>::iterator i = named_groups.begin(); i != named_groups.end(); ++i ){

    NamedGroup* group = (*i).second;
//    if(OPT_VERBOSE){ 
//      std::cout << "Creating volume group " << group->getName() << " of size " << group->getEntities().size() << std::endl;
//    }

    iBase_EntitySetHandle set;
    iGeom_createEntSet(/* 0,*/ &set, &igm_result );
    CHECK_IGEOM( igm_result, "Creating a new entity set " );
    
    const entity_collection_t& group_list = group->getEntities();
    for( entity_collection_t::const_iterator j = group_list.begin(); j != group_list.end(); ++j ){
      iGeom_addEntToSet( *j, set, &igm_result );
      CHECK_IGEOM( igm_result, "Adding entity to material set" );
    }

    std::string name = group->getName();
    /*
    if( name.length() > (unsigned)name_tag_maxlength ){
      name.resize( name_tag_maxlength - 1);
      std::cerr << "Warning: trimmed material name " << group->getName() 
                << " to length " << name_tag_maxlength << std::endl;
    }
    */

    iGeom_setEntSetData( igm, set, name_tag, name.c_str(), &igm_result );
    CHECK_IGEOM( igm_result, "Naming a material group's EntitySet" );
    
  }

}

void MCNP2CAD::GeometryContext::tagCellIDsAsEntNames(){
  
  int igm_result;

  std::string name_tag_id = "NAME";
  int name_tag_maxlength = 64;
  iBase_TagHandle name_tag;

  iGeom_getTagHandle( igm, name_tag_id.c_str(), &name_tag, &igm_result, name_tag_id.length() );
  CHECK_IGEOM( igm_result, "Looking up NAME tag" );
  
  iGeom_getTagSizeBytes( igm, name_tag, &name_tag_maxlength, &igm_result );
  CHECK_IGEOM( igm_result, "Querying NAME tag length" );
//  if( OPT_DEBUG ) std::cout << "Name tag length: " << name_tag_maxlength << " actual id " << name_tag << std::endl;


//  if( OPT_VERBOSE ){ std::cout << "Naming " << named_cells.size() << " volumes." <<  std::endl; }

  for( std::vector< NamedEntity* >::iterator i = named_cells.begin(); i!=named_cells.end(); ++i){
    std::string name = (*i)->getName();
    iBase_EntityHandle entity = (*i)->getHandle();

/*
    if( name.length() > (unsigned)name_tag_maxlength ){
      name.resize( name_tag_maxlength - 1);
      std::cerr << "Warning: trimmed entity name " << (*i)->getName() 
                << " to length " << name_tag_maxlength << std::endl;
    }
    */

    if( entity == NULL ){ std::cerr << "Error: NULL in named_cells" << std::endl; continue; }
    
    iGeom_setData( igm, entity, name_tag, name.c_str(), &igm_result );
    CHECK_IGEOM( igm_result, "Naming an NamedEntity" );
  }

}

bool MCNP2CAD::GeometryContext::defineLatticeNode(  CellCard& cell, iBase_EntityHandle cell_shell, iBase_EntityHandle lattice_shell,
                                          int x, int y, int z, entity_collection_t& accum )
{
  const Lattice& lattice = cell.getLattice();
  int lattice_universe =   cell.getUniverse();

  const FillNode* fn = &(lattice.getFillForNode( x, y, z ));                            
  Transform t = lattice.getTxForNode( x, y, z );        
  int igm_result;
  
  iBase_EntityHandle cell_copy;
  iGeom_copyEnt( cell_shell, &cell_copy, &igm_result );
  CHECK_IGEOM( igm_result, "Copying a lattice cell shell" );
  cell_copy = applyTransform( t, cell_copy );
  
  if( !boundBoxesIntersect( cell_copy, lattice_shell ) ){
    iGeom_deleteEnt( cell_copy, &igm_result);
    CHECK_IGEOM( igm_result, "Deleting a lattice cell shell" );
//    if( OPT_DEBUG ) std::cout << uprefix() << " node failed bbox check" << std::endl;
    return false;
  }

  entity_collection_t node_subcells;
  if( fn->getFillingUniverse() == 0 ){
    // this node of the lattice was assigned universe zero, meaning it's
    // defined to be emtpy. Delete the shell and return true.
    iGeom_deleteEnt( cell_copy, &igm_result );
    CHECK_IGEOM( igm_result, "Deleting a universe-0 lattice cell" );
    return true;
  }
  else if( fn->getFillingUniverse() == lattice_universe ){
    // this node is just a translated copy of the origin element in the lattice
    setVolumeCellID(cell_copy, cell.getIdent());
    if( cell.getMat() != 0 ){ setMaterial( cell_copy, cell.getMat(), cell.getRho() ); }
    if( cell.getImportances().size() ){ setImportances( cell_copy, cell.getImportances()); }
    node_subcells.push_back( cell_copy );
  }
  else{
    // this node has an embedded universe

    iBase_EntityHandle cell_copy_unmoved;
    iGeom_copyEnt( cell_shell, &cell_copy_unmoved, &igm_result );
    CHECK_IGEOM( igm_result, "Re-copying a lattice cell shell" );
    node_subcells = defineUniverse(  fn->getFillingUniverse(), cell_copy_unmoved, (fn->hasTransform() ? &(fn->getTransform()) : NULL ) );
    for( size_t i = 0; i < node_subcells.size(); ++i ){
      node_subcells[i] = applyTransform( t, node_subcells[i] );
    }

    iGeom_deleteEnt( cell_copy, &igm_result );
    CHECK_IGEOM( igm_result, "Deleting lattice cell copy" );

  }

  // bound the node with the enclosing lattice shell
  bool success = false;
  for( size_t i = 0; i < node_subcells.size(); ++i ){
    iBase_EntityHandle lattice_shell_copy;
    iGeom_copyEnt( lattice_shell, &lattice_shell_copy, &igm_result );

    iBase_EntityHandle result;
    if( intersectIfPossible( lattice_shell_copy, node_subcells[i], &result, true ) ){
      updateMaps( node_subcells[i], result );
//      if( OPT_DEBUG ) std::cout << " node defined successfully" << std::endl;
      accum.push_back( result );
      success = true;
    }
    else{ 
      // lattice_shell_copy and node_subcells[i] were deleted by intersectIfPossible(),
      // so there's no need to delete them explicitly
      updateMaps( node_subcells[i], NULL );
//      if( OPT_DEBUG ) std::cout << " node failed intersection" << std::endl;
    }
  }

  return success;
}

MCNP2CAD::MCNP2CAD()// :
//  geom_tag(0), id_tag(0), name_tag(0), category_tag(0), faceting_tol_tag(0), geometry_resabs_tag(0)// Look into this
{
  // set default values
/*  norm_tol = 5;
  faceting_tol = 1e-3;
  len_tol = 0.0;
  */
 // verbose_warnings = false;
//  fatal_on_curves = false;

}

MCNP2CAD::~MCNP2CAD()
{}

std::vector<std::string> MCNP2CAD::get_syntax()
{
  // Define the syntax for the command. Note the syntax is a modified BNF
  // format. Full documentation on the command specification syntax can be
  // found in the documentation.
  std::string syntax =
      "import MCNP2CAD "
      "<string:label='filename',help='<filename>'> "
      "[verbose] [debug] [extra_effort] [skip-mats] "
      "[skip-imps] [skip-nums] [skip_graveyard] "
      "[tol <value:label='specific_tolerance',help='<specific tolerance>'>] ";
      /*
      //input file, igeom_init options?
      "[uwuw_names]";
      /*
      "[faceting_tolerance <value:label='faceting_tolerance',help='<faceting tolerance>'>] "
      "[length_tolerance <value:label='length_tolerance',help='<length tolerance>'>] "
      "[normal_tolerance <value:label='normal_tolerance',help='<normal tolerance>'>] "
      "[verbose] [fatal_on_curves]";
      */

  std::vector<std::string> syntax_list;
  syntax_list.push_back(syntax);

  return syntax_list;
}

std::vector<std::string> MCNP2CAD::get_syntax_help()
{
  std::vector<std::string> help;
  return help;
}

std::vector<std::string> MCNP2CAD::get_help()
{
  std::vector<std::string> help;
  return help;
}

bool MCNP2CAD::execute(CubitCommandData &data)
{

  std::string mcnp2cad_version(bool full = true);

  struct program_option_struct Gopt;

  Gopt.verbose = Gopt.debug = false;
  Gopt.infinite_lattice_extra_effort = false;
  Gopt.tag_materials = true;
  Gopt.tag_importances = true;
  Gopt.tag_cell_IDs = true;
  Gopt.make_graveyard = true;
//  Gopt.input_file = "";
//  Gopt.igeom_init_options = "";
  Gopt.override_tolerance = false;
  Gopt.uwuw_names = false;

  bool DiFlag = false, DoFlag = false;

  ProgOptions po("mcnp2cad " + mcnp2cad_version(false) +  ": An MCNP geometry to CAD file converter");
  po.setVersion( mcnp2cad_version() );

  po.addOpt<void>("extra-effort,e","Use extra effort to get infinite lattices right (may be slow)", 
                  &Gopt.infinite_lattice_extra_effort );
  po.addOpt<void>("verbose,v", "Verbose output", &Gopt.verbose, po.store_true );
  po.addOpt<void>("debug,D", "Debugging (very verbose) output", &Gopt.debug, po.store_true );
  po.addOpt<void>("Di", "Debug output for MCNP parsing phase only", &DiFlag, po.store_true );
  po.addOpt<void>("Do","Debug output for iGeom output phase only", &DoFlag, po.store_true );

  po.addOptionHelpHeading( "Options controlling CAD output:" );
  po.addOpt<double>("tol,t", "Specify a tolerance for merging surfaces", &Gopt.specific_tolerance );
  po.addOpt<void>("skip-mats,M", "Do not tag materials using group names", 
                  &Gopt.tag_materials, po.store_false );
  po.addOpt<void>("skip-imps,P", "Do not tag cell importances using group names",
                  &Gopt.tag_importances, po.store_false );
  po.addOpt<void>("skip-nums,N", "Do not tag cell numbers using body names",
                  &Gopt.tag_cell_IDs, po.store_false );
  po.addOpt<void>("skip-graveyard,G", "Do not bound the geometry with a `graveyard' bounding box",
                  &Gopt.make_graveyard, po.store_false );
  po.addOpt<void>("uwuw-names,U", "Use a UWUW compatible name scheme for material groups,"
                                   "i.e. 'mat:mX/rho:Y' where X is material number is Y is density",
                  &Gopt.uwuw_names, po.store_true );

  po.addRequiredArg( "input_file", "Path to MCNP geometry input file", &Gopt.input_file );


  std::string filename;
  data.get_string("filename",filename);

  std::ifstream input( filename.c_str(), std::ios::in );
  InputDeck& deck = InputDeck::build(input);
  std::string output = "filename is " + filename + "\n";
  PRINT_INFO( output.c_str() );

  iGeom_Instance igm;
  int igm_result;

  iGeom_newGeom( Gopt.igeom_init_options.c_str(), &igm, &igm_result, Gopt.igeom_init_options.length() );
  CHECK_IGEOM( igm_result, "Initializing iGeom" );

  GeometryContext context( igm, deck );
  context.createGeometry();
  
  /*
  mdbImpl = new moab::Core();
  myGeomTool = new moab::GeomTopoTool(mdbImpl);
  message.str("");

  bool result = true;
  moab::ErrorCode rval;

  // Create entity sets for all geometric entities
  refentity_handle_map entmap[5];

  rval = create_tags();
  CHK_MB_ERR_RET("Error initializing DAGMC export: ",rval);

  rval = parse_options(data);
  CHK_MB_ERR_RET("Error parsing options: ",rval);
  
  // Always tag with the faceting_tol and geometry absolute resolution
  // If file_set is defined, use that, otherwise (file_set == NULL) tag the interface
  moab::EntityHandle set = 0;
  rval = mdbImpl->tag_set_data(faceting_tol_tag, &set, 1, &faceting_tol);
  CHK_MB_ERR_RET("Error tagging geometry with faceting tolerance: ", rval)

  rval = mdbImpl->tag_set_data(geometry_resabs_tag, &set, 1, &GEOMETRY_RESABS);
  CHK_MB_ERR_RET("Error tagging geometry with absolute geometry resolution: ", rval)

  rval = create_entity_sets(entmap);
  CHK_MB_ERR_RET("Error creating entity sets: ",rval);

  rval = create_topology(entmap);
  CHK_MB_ERR_RET("Error creating topology: ",rval);
  
  rval = store_surface_senses(entmap[2], entmap[3]);
  CHK_MB_ERR_RET("Error storing surface senses: ",rval);
  
  rval = store_curve_senses(entmap[1], entmap[2]);
  CHK_MB_ERR_RET("Error storing curve senses: ",rval);
    
  rval = store_groups(entmap);
  CHK_MB_ERR_RET("Error storing groups: ",rval);
  
  entmap[3].clear();
  entmap[4].clear();
  
  rval = create_vertices(entmap[0]);
  CHK_MB_ERR_RET("Error creating vertices: ",rval);
  
  rval = create_curve_facets(entmap[1], entmap[0]);
  CHK_MB_ERR_RET("Error faceting curves: ",rval);

  rval = create_surface_facets(entmap[2], entmap[0]);
  CHK_MB_ERR_RET("Error faceting surfaces: ",rval);

*/
  /*
  rval = mdbImpl->write_file(filename.c_str());
  CHK_MB_ERR_RET("Error writing file: ",rval);

  teardown();
  
  return result;
  */
  return true;
}


std::string mcnp2cad_version( bool full ){
  std::stringstream str;
  str << (full ? "mcnp2cad version " : "")
      << MCNP2CAD_VERSION_MAJOR << "." 
      << MCNP2CAD_VERSION_MINOR << "." 
      << MCNP2CAD_VERSION_REV;
  if(full)
      str << "\nCompiled on " << __DATE__ << " at " << __TIME__ ;
  return str.str();
}
//}

