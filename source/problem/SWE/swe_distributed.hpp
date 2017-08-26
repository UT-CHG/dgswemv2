#ifndef SWE_DISTRIBUTED_HPP
#define SWE_DISTRIBUTED_HPP

#include <vector>

namespace SWE {
class Distributed {
  private:
    std::vector<double>& send_buffer;
    std::vector<double>& receive_buffer;

  public:
    Distributed(std::vector<double>& send_buf, std::vector<double>& recv_buf)
        : send_buffer(send_buf), receive_buffer(recv_buf) {}

    inline std::vector<double>::iterator GetSendBuffIterator(uint index) {
        return send_buffer.begin() + index;
    }

    inline std::vector<double>::const_reverse_iterator GetReceiveBufferRIterator(uint indx) {
        return receive_buffer.crbegin() + indx;
    }
};
}
#endif