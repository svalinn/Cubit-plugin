#ifndef CATAG_HPP
#define CATAG_HPP

#include "iGeom.h"
#include "RefGroup.hpp"
/**\file CATag.hpp
 *
 *\author Tim Tautges
 *\author Jason Kraftcheck
 *
 * Original file from SNL TSTT repository was named CATSTT.
 *
 * Renamed CATag and added to ANL ITAPS repository by J.Kraftcheck,
 * 2007-6-15
 */

class CATag;

class CGMTagManager
{
public:
#ifdef ITAPS_SHIM
//  iGeom_vtable *vtable;
#endif

  friend class CATag;

//  If the destructor is there, the plugin isn't recognized.
//  ~CGMTagManager();

  struct TagInfo
  {
    int tagLength;
    std::string tagName;
    int tagType;
    char *defaultValue;
    bool isActive;
  };

  static CubitAttrib* CATag_creator(RefEntity* entity, CubitSimpleAttrib &p_csa);

  iBase_ErrorType createTag (/*in*/ const char *tag_name,
                             /*in*/ const int tag_size,
                             /*in*/ const int tag_type,
                             /*in*/ char* default_value,
                             /*out*/ long *tag_handle);

  int getTagSize (/*in*/ const long tag_handle);

  long getTagHandle (/*in*/ const char *tag_name);

  iBase_ErrorType setArrData (/*in*/ RefEntity* const* entity_handles,
                                     const int entity_handles_size,
                              /*in*/ const long tag_handle,
                              /*in*/ const char *tag_values);

  static inline CGMTagManager& instance()
  {
    static CGMTagManager static_instance;
    return static_instance;
  }

private:
  CGMTagManager();

  int CATag_att_type;

  std::vector<TagInfo> tagInfo;

  static TagInfo* const presetTagInfo;
  static const int numPresetTag;

  std::map<std::string, long> tagNameMap;
  static const char *CATag_NAME;
  static const char *CATag_NAME_INTERNAL;

  RefGroup *interfaceGroup;

  iBase_ErrorType setPresetTagData(RefEntity *entity, const long tag_num,
                                   const char *tag_value, const int tag_size);

  CATag *get_catag(RefEntity *ent,
                   const bool create_if_missing = false);

  RefGroup *interface_group(const bool create_if_missing = true);

};

class CATag: public CubitAttrib
{
private:
  friend class CGMTagManager;

  std::map<int, void*> tagData;

  CGMTagManager *myManager;

  CATag(CGMTagManager *manager, RefEntity *owner, CubitSimpleAttrib *csa_ptr);
    //- create a CATag from a simple attribute

public:

  CubitStatus actuate() {return CUBIT_SUCCESS;}

//  the destructor will make the plugin not be recognized.
//  virtual ~CATag();

  CubitStatus update();

  CubitStatus reset();

  CubitSimpleAttrib cubit_simple_attrib();

  int int_attrib_type() {return myManager->CATag_att_type;}

  void add_csa_data(CubitSimpleAttrib *csa_ptr);

//  void print();

  iBase_ErrorType set_tag_data(long tag_num, const void *tag_data,
                               const bool can_shallow_copy = false);
};

static void tokenize(const std::string& str,
                     std::vector<std::string>& tokens);

// Expect option of the form "NAME=VALUE".
// If NAME portion matches, pass back VALUE and return true.
// Otherwise, leave 'value' unchanged and return false.
static bool match_option(const std::string& opt,
                         const char* name,
                         std::string& value);

CubitStatus CubitCompat_export_solid_model(DLIList<RefEntity*>& ref_entity_list,
                                           const char* filename,
                                           const char * filetype,
                                           int &num_ents_exported,
                                           const CubitString &cubit_version,
                                           const char* logfile_name = NULL);

static CubitStatus iGeom_bounding_box(RefEntity* entity,
                                      CubitVector& minc,
                                      CubitVector& maxc);

static void copy_ibase_type(int ibase_type,
                            const DLIList<CubitEntity*>& list,
                            iBase_EntityHandle** entity_handles,
                            int* entith_handles_alloc,
                            int* entity_handles_size,
                            int* err);

static void append_all_ibase_type(int ibase_type,
                                  DLIList<RefEntity*>& target_list,
                                  int* err);

static int count_ibase_type(int ibase_type,
                            const DLIList<CubitEntity*>& list,
                            int* err);

template<typename SKIP_TYPE> static
int append_not_type(const DLIList<CubitEntity*>& source_list,
                    iBase_EntityHandle* array,
                    int array_size);

template <typename TARGET_TYPE> static
int append_type(const DLIList<CubitEntity*>& source_list,
                iBase_EntityHandle* array,
                int array_size);

template <typename T> static
int count_type(const DLIList<CubitEntity*>& list);
#endif //CATAG_HPP
