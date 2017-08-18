#ifndef HPX_COMMUNICATOR_HPP
#define HPX_COMMUNICATOR_HPP

HPX_REGISTER_CHANNEL(uint);

template <typename ProblemType>
class HPXCommunicator {
  private:
    hpx::lcos::channel<uint> outgoing;
    hpx::lcos::channel<uint> incoming;

  public:
    HPXCommunicator() = default;
    HPXCommunicator(std::string, std::string);

    void Send(uint, double);
    hpx::future<uint> Receive(double);
};

template <typename ProblemType>
HPXCommunicator<ProblemType>::HPXCommunicator(std::string in_location, std::string out_location) {
    outgoing = hpx::lcos::channel<uint>(hpx::find_here());
    incoming = hpx::lcos::channel<uint>();

    std::string outgoing_channel_string = "channel_from_" + in_location + "_to_" + out_location;
    std::string incoming_channel_string = "channel_from_" + out_location + "_to_" + in_location;

    hpx::future<void> set_out = outgoing.register_as(outgoing_channel_string);
    incoming.connect_to(incoming_channel_string);

    set_out.get();
}

template <typename ProblemType>
void HPXCommunicator<ProblemType>::Send(uint message, double timestamp) {
    outgoing.set(message, timestamp);
}

template <typename ProblemType>
hpx::future<uint> HPXCommunicator<ProblemType>::Receive(double timestamp) {
    return incoming.get(timestamp);
}

#endif