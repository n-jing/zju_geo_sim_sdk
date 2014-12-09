#ifndef PROPERTY_SXX_H
#define PROPERTY_SXX_H
#include <zjucad/matrix/matrix.h>
#include <cassert>
#include <vector>
#include <set>
#include <string>

namespace sxx
{
  class Property
  {
  public:
    Property():property_name_("property")
    {
    }

    Property(const std::string &str):property_name_(str)
    {
    }

    const std::string &get_pro_name() const
    {
      return property_name_;
    }

    virtual int cal_property(const zjucad::matrix::matrix<double> &weight,
                            const zjucad::matrix::matrix<double> &value)
    {
      return 0;
    }

  private:
    std::string property_name_;

  };

  class Frame_property:public Property
  {
  public:
    Frame_property():Property("frame_property"),frame_value_(3, 1)
    {
    }

    Frame_property(const zjucad::matrix::matrix<double> &value):Property("frame_property")
    {
      assert(value.size(1) == 3 && value.size(2) == 1);
      frame_value_ = value;
    }

    //  @param weight is 3 * 1
    //  @param value is 3 * 3
    virtual int cal_property(const zjucad::matrix::matrix<double> &weight,
                            const zjucad::matrix::matrix<double> &value)
    {
      assert(get_pro_name() == "frame_property");
      assert(weight.size(1) == 3 && weight.size(2) == 1);
      assert(value.size(1) == 3 && value.size(2) == 3);
      frame_value_ = value * weight;
    }

    const zjucad::matrix::matrix<double> &get_frame() const
    {
      return frame_value_;
    }

    zjucad::matrix::matrix<double> &get_frame()
    {
      return frame_value_;
    }

    void set_frame(const zjucad::matrix::matrix<double> &value)
    {
      assert(value.size(1) == 3 && value.size(2) == 1);
      frame_value_ = value;
    }

  private:
    //frame_value_ is 3*1
    zjucad::matrix::matrix<double> frame_value_;
  };

  class UV_property:public Property
  {
  public:
    UV_property():Property("uv_property"),uv_value_(2, 1)
    {
    }

    UV_property(const zjucad::matrix::matrix<double> &value):Property("uv_property")
    {
      assert(value.size(1) == 2 && value.size(2) == 1);
      uv_value_ = value;
    }

    void set_uv_value(const zjucad::matrix::matrix<double> &value)
    {
      uv_value_ = value;
    }

    //  @param weight is 3 * 1
    //  @param value is 2 * 3
    virtual int cal_property(const zjucad::matrix::matrix<double> &weight,
                            const zjucad::matrix::matrix<double> &value)
    {
      assert(get_pro_name() == "uv_property");
      assert(weight.size(1) == 3 && weight.size(2) == 1);
      assert(value.size(1) == 2 && value.size(2) == 3);
      uv_value_ = value * weight;
      return 0;
    }

    const zjucad::matrix::matrix<double> &get_uv() const
    {
      return uv_value_;
    }

  private:
    //  uv_value_ is 2 * 1
    zjucad::matrix::matrix<double> uv_value_;
  };
}

#endif // PROPERTY_SXX_H
