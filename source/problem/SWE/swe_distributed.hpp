#ifndef SWE_DISTRIBUTED_HPP
#define SWE_DISTRIBUTED_HPP

#include <vector>

namespace SWE {
class Distributed {
  private:
    std::vector<double>& send_buffer;
    std::vector<double>& receive_buffer;

    uint send_ze_indx;
    uint send_qx_indx;
    uint send_qy_indx;

    uint recv_ze_indx;
    uint recv_qx_indx;
    uint recv_qy_indx;

  public:
    Distributed(std::vector<double>& send_buf,
                std::vector<double>& recv_buf,
                uint send_ze_indx,
                uint send_qx_indx,
                uint send_qy_indx,
                uint recv_ze_indx,
                uint recv_qx_indx,
                uint recv_qy_indx)
        : send_buffer(send_buf),
          receive_buffer(recv_buf),
          send_ze_indx(send_ze_indx),
          send_qx_indx(send_qx_indx),
          send_qy_indx(send_qy_indx),
          recv_ze_indx(recv_ze_indx),
          recv_qx_indx(recv_qx_indx),
          recv_qy_indx(recv_qy_indx) {}

    inline std::vector<double>::iterator GetSendZeIterator() { return send_buffer.begin() + send_ze_indx; }

    inline std::vector<double>::iterator GetSendQxIterator() { return send_buffer.begin() + send_qx_indx; }

    inline std::vector<double>::iterator GetSendQyIterator() { return send_buffer.begin() + send_qy_indx; }

    inline std::vector<double>::const_reverse_iterator GetReceiveZeRIterator() {
        return receive_buffer.crbegin() + recv_ze_indx;
    }

    inline std::vector<double>::const_reverse_iterator GetReceiveQxRIterator() {
        return receive_buffer.crbegin() + recv_qx_indx;
    }

    inline std::vector<double>::const_reverse_iterator GetReceiveQyRIterator() {
        return receive_buffer.crbegin() + recv_qy_indx;
    }
};
}
#endif