#ifndef GMESH_MESH_ITEM_STATUS_H_
#define GMESH_MESH_ITEM_STATUS_H_

namespace gmesh{

enum status_bits {
  DELETED      = 1,   ///< item has been deleted
  UNUSED       = 256 ///< unused
};

class status_info {
public:
  typedef unsigned int value_type;
  status_info() : status_(0) {}

  bool is_deleted() const { return is_bit_set(DELETED); }
  void set_deleted(bool b) { change_bit(DELETED, b); }
  
  unsigned int get_bits() const { return status_; }
  void set_bits(unsigned int bits) { status_ = bits; }

  bool is_bit_set(unsigned int s) const { return (status_ & s) > 0; }
  void set_bit(unsigned int s) { status_ |= s; }
  void unset_bit(unsigned int s) { status_ &= ~s; }
  void change_bit(unsigned int s, bool b) {
    if(b) status_ |= s;
    else status_ &= ~s;
  }
    
private:
  value_type status_;
};
}

#endif
