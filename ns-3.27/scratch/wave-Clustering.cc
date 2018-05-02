/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2014 Dalian University of Technology
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR average_throughput PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received average_throughput copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Author: Junling Bu <linlinjavaer@gmail.com>
 */
#include "ns3/node.h"
#include "ns3/packet.h"
#include "ns3/simulator.h"
#include "ns3/log.h"
#include "ns3/node-container.h"
#include "ns3/net-device-container.h"
#include "ns3/yans-wifi-helper.h"
#include "ns3/mobility-helper.h"
#include "ns3/rng-seed-manager.h"
#include "ns3/command-line.h"
#include "ns3/mobility-model.h"
#include <iostream>
#include <algorithm>
#include "ns3/string.h"
#include "ns3/double.h"
#include "ns3/header.h"
#include "ns3/object-map.h"
#include "ns3/regular-wifi-mac.h"
#include "ns3/constant-velocity-mobility-model.h"
#include "ns3/wave-net-device.h"
#include "ns3/wave-mac-helper.h"
#include "ns3/wave-helper.h"
#include "ns3/netanim-module.h"
#include<cstdlib>
#include<ctime>
#include <sys/time.h>
#include <stdlib.h>
#include <unistd.h>

using namespace ns3;

NS_LOG_COMPONENT_DEFINE ("WaveMultipleChannel");

class StatsTag : public Tag
{
public:
  StatsTag (void)
    : m_packetId (0),
      m_sendTime (Seconds (0)),
      m_packetType(0),
	  m_nodeId(0)
  {
  }
  StatsTag (uint32_t packetId, Time sendTime,uint32_t pactetType,uint32_t nodeid)
    : m_packetId (packetId),
      m_sendTime (sendTime),
	  m_packetType(pactetType),
	  m_nodeId(nodeid)
  {
  }
  virtual ~StatsTag (void)
  {
  }

  uint32_t GetPacketId (void)
  {
    return m_packetId;
  }
  Time GetSendTime (void)
  {
    return m_sendTime;
  }

  static TypeId GetTypeId (void);
  virtual TypeId GetInstanceTypeId (void) const;
  virtual uint32_t GetSerializedSize (void) const;
  virtual void Serialize (TagBuffer i) const;
  virtual void Deserialize (TagBuffer i);
  virtual void Print (std::ostream &os) const;

	uint32_t getPacketType() const {
		return m_packetType;
	}

	void setPackettType(uint32_t packetType) {
		m_packetType = packetType;
	}

	uint32_t getNodeId() const {
		return m_nodeId;
	}

	void setNodeId(uint32_t nodeId) {
		this->m_nodeId = nodeId;
	}

private:
  uint32_t m_packetId;
  Time m_sendTime;
  uint32_t m_packetType;
  uint32_t m_nodeId;
};
TypeId
StatsTag::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::StatsTag")
    .SetParent<Tag> ()
    .AddConstructor<StatsTag> ()
  ;
  return tid;
}
TypeId
StatsTag::GetInstanceTypeId (void) const
{
  return GetTypeId ();
}
uint32_t
StatsTag::GetSerializedSize (void) const
{
  return 3*sizeof (uint32_t) + sizeof (uint64_t);
}
void
StatsTag::Serialize (TagBuffer i) const
{
  i.WriteU32 (m_packetId);
  i.WriteU64 (m_sendTime.GetMicroSeconds ());
  i.WriteU32 (m_packetType);
  i.WriteU32(m_nodeId);
}
void
StatsTag::Deserialize (TagBuffer i)
{
  m_packetId = i.ReadU32 ();
  m_sendTime = MicroSeconds (i.ReadU64 ());
  m_packetType = i.ReadU32 ();
  m_nodeId=  i.ReadU32 ();
}
void
StatsTag::Print (std::ostream &os) const
{
  os << "packet=" << m_packetId << " sendTime=" << m_sendTime<< "  packetType="<<m_packetType;
}


class NodeInfoHeader :public Header{
public:
	NodeInfoHeader();
	virtual ~NodeInfoHeader ();

	static TypeId GetTypeId (void);
	virtual TypeId GetInstanceTypeId (void) const;
	virtual void Print (std::ostream &os) const;
	virtual void Serialize (Buffer::Iterator start) const;
	virtual uint32_t Deserialize (Buffer::Iterator start);
	virtual uint32_t GetSerializedSize (void) const;

	double getSpeed() const {
		return m_speed;
	}

	void setSpeed(double speed) {
		m_speed = speed;
	}

	const Vector& getPosition() const {
		return m_position;
	}

	void setPosition(const Vector& position) {
		m_position = position;
	}

	double getRiskFactor() const {
		return m_riskFactor;
	}

	void setRiskFactor(double riskFactor) {
		m_riskFactor = riskFactor;
	}

private:
	Vector m_position;
	double m_riskFactor;
	double m_speed;
};
NodeInfoHeader::NodeInfoHeader(){}
NodeInfoHeader::~NodeInfoHeader(){}
TypeId
NodeInfoHeader::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::NodeInfoHeader")
    .SetParent<Header> ()
    .AddConstructor<NodeInfoHeader> ()
  ;
  return tid;
}
TypeId
NodeInfoHeader::GetInstanceTypeId (void) const
{
  return GetTypeId ();
}

void
NodeInfoHeader::Print (std::ostream &os) const
{
  os << "position=" << m_position<<" risFactor="<<m_riskFactor<<" speed="<<m_speed<<std::endl;;
}
uint32_t
NodeInfoHeader::GetSerializedSize (void) const
{
  // we reserve 2 bytes for our header.
  return 5*sizeof(double);
}
void
NodeInfoHeader::Serialize (Buffer::Iterator start) const
{
  // we can serialize two bytes at the start of the buffer.
  // we write them in network byte order.
  //start.WriteHtonU16 (m_data);
	double x = m_position.x;
	start.WriteDouble(x);
	double y =m_position.y;
	start.WriteDouble(y);
	double z =m_position.z;
	start.WriteDouble(z);
	start.WriteDouble(m_riskFactor);
	start.WriteDouble(m_speed);

}
uint32_t
NodeInfoHeader::Deserialize (Buffer::Iterator start)
{
	double x = start.ReadDouble();
	double y = start.ReadDouble();
	double z = start.ReadDouble();
	m_position = Vector(x,y,z);
	m_riskFactor = start.ReadDouble();
	m_speed = start.ReadDouble();
	return 5*sizeof(double);
}

struct Cdq{
	uint32_t send_id;
	double comend;
};


class CdqHeader :public Header{
public:
	CdqHeader();
	virtual ~CdqHeader ();
	static TypeId GetTypeId (void);
	virtual TypeId GetInstanceTypeId (void) const;
	virtual void Print (std::ostream &os) const;
	virtual void Serialize (Buffer::Iterator start) const;
	virtual uint32_t Deserialize (Buffer::Iterator start);
	virtual uint32_t GetSerializedSize (void) const;

	const std::vector<struct Cdq> getCdqVector() const {
		return cdq_vector;
	}

	void setCdqVector(const std::vector<struct Cdq> cdqVector) {
		cdq_vector = cdqVector;
	}

	uint32_t getSize() const {
		return size;
	}

	void setSize(uint32_t size) {
		this->size = size;
	}

	uint32_t getChannel() const {
		return channel;
	}

	void setChannel(uint32_t channel) {
		this->channel = channel;
	}


private:
	uint32_t channel;
	uint32_t size;
	std::vector<struct Cdq> cdq_vector;
};
CdqHeader::CdqHeader(){}
CdqHeader::~CdqHeader(){}
TypeId
CdqHeader::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::CdqHeader")
    .SetParent<Header> ()
    .AddConstructor<CdqHeader> ()
  ;
  return tid;
}
TypeId
CdqHeader::GetInstanceTypeId (void) const
{
  return GetTypeId ();
}

void
CdqHeader::Print (std::ostream &os) const
{
 // os << "position=" << m_position<<" risFactor="<<m_riskFactor<<" speed="<<m_speed<<std::endl;;
}
uint32_t
CdqHeader::GetSerializedSize (void) const
{
	  return size*(sizeof(double)+sizeof(uint32_t)) +2*sizeof(uint32_t);
}
void
CdqHeader::Serialize (Buffer::Iterator start) const
{
	start.WriteU32(channel);
	start.WriteU32(size);
	struct Cdq cdq;
	for(int i=0;i<size;i++){
		cdq=cdq_vector[i];
		start.WriteU32(cdq.send_id);
		start.WriteDouble(cdq.comend);
	}
}
uint32_t
CdqHeader::Deserialize (Buffer::Iterator start)
{
	channel = start.ReadU32();
	size = start.ReadU32();
	for(int i=0;i<size;i++){
		struct Cdq cdq;
		cdq.send_id = start.ReadU32();
		cdq.comend = start.ReadDouble();
		cdq_vector.push_back(cdq);
	}
  return size*(sizeof(double)+sizeof(uint32_t)) +2*sizeof(uint32_t);
}

// it is useless here we can use any value we like
// because we only allow to send and receive one type packets
#define Packet_Number 0x8888

// the transmission range of a device should be decided carefully,
// since it will affect the packet delivery ratio
// we suggest users use wave-transmission-range.cc to get this value
#define Device_Transmission_Range 50
#define Max_timeslot 10  //the max number time slot
#define Car_number 20  //the mumber of cars

//PacketType
#define RmPacket 0
#define CdqPacket 1
#define AmPacket 2
#define SchPacket 3




//the of sumilator start


/**
 * This simulation is mainly used for wns3 to evaluate the performance of IEEE 1609.4
 * with current implementation.
 *
 * This simulation contains many parameters we can configure such as the number of nodes.
 * In our paper,
 * First we study the impact of the standard on safety packets.
 * We use small packet size (default 200) and broadcast transmission mechanism to model safety packets.
 * we consider four cases of the standard: A-not use, B-normal usage, C-best usage, D-worst usage.
 * And we change the number of nodes to study the communication performance.
 * Second we study the impact of the standard on non-safety packets.
 * We use large packet size (default 1500) and unicast transmission mechanism to model non-safety packets.
 * we also consider four cases of the standard: A-not use, B-normal usage, C-best usage, D-worst usage.
 * And we also change the number of nodes to study the communication performance. *
 */
class MultipleChannelsExperiment
{
public:
  MultipleChannelsExperiment (void);

  bool Configure (int argc, char **argv);
  void Usage (void);
  //get current time ，unit: s
  double getCurrentTime();
  void Run (void);

private:

  void CreateNodes (void);
  void CreateWaveDevice (void);
  void SetupMobility (void);

  void InstallApplicationA (void);
  void SendRmPackets(Time time, uint32_t channelNumber);
  void SendCdqPackets(uint32_t channelNumber);
  void SendAmPackets(Ptr<WaveNetDevice> sender, uint32_t channelNumber,bool falg);
  void SendSchPackets(bool flag);
  void SendSch(Time time);

  void SendCdq(Ptr<WaveNetDevice> sender);
  void SendRm(Ptr<WaveNetDevice> sender,uint32_t channelNumber,bool flag);
  void Send (Ptr<WaveNetDevice> sender, uint32_t channelNumber);
  bool Receive (Ptr<NetDevice> dev, Ptr<const Packet> pkt, uint16_t mode, const Address &sender);

  void RmPacketEvent(Ptr<NetDevice> dev, Ptr<const Packet> pkt, uint32_t send_id);
  void CdqPacketEvent(Ptr<NetDevice> dev, Ptr<const Packet> pkt, const Address &sender);
  void AmPacketEvent(Ptr<NetDevice> dev, Ptr<const Packet> pkt, const Address &sender);
  void SchPacketEvent(Ptr<NetDevice> dev, Ptr<const Packet> pkt, const Address &sender);

  void InitStats (void);
  void Destroy (void);
  void Stats (uint32_t randomNumber);
  void StartSch();
  void StopSch();

  //void StatQueuedPackets (void);
 // uint32_t GetQueuedSize (Ptr<WaveEdcaTxopN> edca);
  uint32_t GetChannel();


  NodeContainer nodes;
  NetDeviceContainer devices;

  uint32_t nodesNum;          // current nodes in 4-lines and 1km
  uint32_t freq;              // the frequency for sending packets
  uint32_t simulationTime;    // simulation time
  uint32_t size;              // the size of packet
  uint32_t rm_size;           //the size of rm packet
  //uint32_t cdq_size;           //the size of cdq packet
  //uint32_t am_size;           //the size of am packet

  double time_start=0.0;

  uint32_t Channel[6] = {SCH1,SCH2,SCH3,SCH4,SCH5,SCH6};
  uint32_t channel_current=0;
  uint32_t start_chanel[Car_number];

  Ptr<UniformRandomVariable> rng;

  // used for identifying received packet
  uint32_t sequence;

  // every received packets will be used for stat the duplicate packet
  // and duplicate packet case will only happen in unicast
  std::vector <uint32_t> receivedPackets;
  // we will check whether the packet is received by the neighbors
  // in the transmission range
  std::map<uint32_t, std::vector<uint32_t> *> broadcastPackets;
  //std::map<uint32_t, std::vector<struct CdqInfo>*> RMs;
  std::map<uint32_t, std::vector<struct Cdq> *> CDQ;
  std::map<uint32_t, std::vector<struct Cdq> *> CDQ1;
  std::map<uint32_t, std::vector<uint32_t> *> CLUSTER;
  uint32_t channel_table[Car_number]={172};
  Address nodes_cluster_address[Car_number];
  uint32_t nodes_cluster[Car_number];
  bool isCluster[Car_number];
  bool isSuccefulCluster[Car_number];
  double distant_from_cluster[Car_number]={0.0};
  double risk_factor[Car_number]={0.0};
  uint32_t startsch_chanel[Car_number]={172};

/*
  class SendStat
  {
  public:
    uint32_t inCchi;
    uint32_t inCgi;
    uint32_t inSchi;
    uint32_t inSgi;
    SendStat (): inCchi (0), inCgi (0), inSchi (0), inSgi (0){}
  };
  */
  uint32_t success_cluster;
  uint32_t unsuccess_cluster;
  uint32_t sends;
  uint32_t receives;
  uint64_t delaySum;        // us
  uint64_t receiveSum;

  uint32_t run;
  double pdr,delay,system_throughput,average_throughput;
};

MultipleChannelsExperiment::MultipleChannelsExperiment (void)
  : nodesNum (Car_number),           // 20 vehicles in 1km * 4-line
    freq (10),               // 10Hz, 100ms send one non-safety packet
    simulationTime (10),     // make it run 100s
    size (1500),             // packet size
	rm_size(200),
	//cdq_size(1500),
	//am_size(200),
    sequence (0),
    run (1),
    pdr (0),
    delay (0),
    system_throughput (0),
    average_throughput (0)
{
}

bool
MultipleChannelsExperiment::Configure (int argc, char **argv)
{
  CommandLine cmd;
  cmd.AddValue ("nodes", "Number of nodes.", nodesNum);
  cmd.AddValue ("time", "Simulation time, s.", simulationTime);
  cmd.AddValue ("size", "Size of safety packet, bytes.", size);
  cmd.AddValue ("frequency", "Frequency of sending safety packets, Hz.", freq);
  cmd.AddValue ("run", "to make the result more credible, make simulation run many times with different random number).", run);

  cmd.Parse (argc, argv);
  return true;
}

void
MultipleChannelsExperiment::Usage (void)
{
  std::cout << "usage:"
		    << "./waf --run=\"wave-multiple-channel --nodes=20 --time=100 --size=200 "
		    		"--frequency=10 --broadcast=false \""
            << std::endl;
}

double MultipleChannelsExperiment::getCurrentTime() {
	struct timeval tv;
	gettimeofday(&tv, NULL);    //该函数在sys/time.h头文件中
	return tv.tv_sec + tv.tv_usec / 1000000.0;
}


void
MultipleChannelsExperiment::CreateNodes ()
{
  NS_LOG_FUNCTION (this);
  NS_ASSERT (nodesNum != 0);
  nodes = NodeContainer ();
  nodes.Create (nodesNum);
  NS_LOG_DEBUG ("create " << nodesNum << " nodes");
}

void
MultipleChannelsExperiment::CreateWaveDevice (void)
{
  NS_LOG_FUNCTION (this);
  // the transmission range is about 250m refer to "wave-transmission-range.cc"
  YansWifiChannelHelper wifiChannel;
  wifiChannel.SetPropagationDelay ("ns3::ConstantSpeedPropagationDelayModel");
  wifiChannel.AddPropagationLoss ("ns3::TwoRayGroundPropagationLossModel",
 			  	  	  	  	  	  "Frequency", DoubleValue (5.89e9),
 			  	  	  	  	      "HeightAboveZ", DoubleValue (0.5));
  YansWifiPhyHelper wifiPhy =  YansWifiPhyHelper::Default ();
  wifiPhy.Set("TxPowerStart",  DoubleValue (16));
  wifiPhy.Set("TxPowerEnd",  DoubleValue (16));
  wifiPhy.SetChannel (wifiChannel.Create ());
  QosWaveMacHelper waveMac = QosWaveMacHelper::Default ();
  WaveHelper waveHelper = WaveHelper::Default ();
  devices = waveHelper.Install (wifiPhy, waveMac, nodes);

  // Enable WAVE logs.
  // WaveHelper::LogComponentEnable();
  // or waveHelper.LogComponentEnable();

  for (uint32_t i = 0; i != devices.GetN (); ++i)
      {
      devices.Get (i)->SetReceiveCallback (MakeCallback (&MultipleChannelsExperiment::Receive,this));
      }
}

void
MultipleChannelsExperiment::SetupMobility ()
{
  // Create a simple highway model based on GridPositionAllocator and ConstantPositionMobilityModel
  MobilityHelper mobility;
  // every node is in 1km * 4line. It is 3 meter between neighboring lines.
  mobility.SetPositionAllocator ("ns3::GridPositionAllocator",
		  	  	  	  	  	  	 "MinX", DoubleValue (0.0),
		  	  	  	  	  	  	 "MinY", DoubleValue (0.0),
		  	  	  	  	  	  	 "DeltaX", DoubleValue (1000.0 / (nodesNum / 4.0)),
		  	  	  	  	  	  	 "DeltaY", DoubleValue (5),
		  	  	  	  	  	  	 "GridWidth", UintegerValue (nodesNum / 4),
		  	  	  	  	  	  	 "LayoutType", StringValue ("RowFirst"));
  // here we use constant velocity mobility model for two reasons
  // (a) since the mobility model of ns3 is mainly useful and enough for MANET, we can not ensure
  // whether they are suitable for VANET. Some special mobility models are required for VANET simulation.
  // (b) the mobility characteristic can cause packets lost which we do not want here. Because we want
  // get the impact of the standard on communication performance.
  mobility.SetMobilityModel ("ns3::ConstantVelocityMobilityModel");
  //mobility.SetMobilityModel ("ns3::ConstantPositionMobilityModel");
  mobility.Install (nodes);
  srand((unsigned)time(NULL));
  for (NodeContainer::Iterator i = nodes.Begin (); i != nodes.End (); ++i)
    {
	  	  Ptr<Node> node = (*i);
	  	  node->GetObject<ConstantVelocityMobilityModel> ()->SetVelocity (Vector (15+rand()%Car_number, 0, 0));
	  	  Ptr<MobilityModel> model = node->GetObject<MobilityModel> ();
	  	  Vector pos = model->GetPosition ();
	  	  NS_LOG_DEBUG ( "node :" << node->GetId() << " position: " << pos);
    }
}
void
MultipleChannelsExperiment::RmPacketEvent(Ptr<NetDevice> dev, Ptr<const Packet> pkt, uint32_t send_id){
	StatsTag tag;
	bool result;
	result = pkt->FindFirstMatchingByteTag (tag);
	if (!result){
		NS_FATAL_ERROR ("the packet here shall have a stats tag");
	}
		 // std::cout<<"packet :"<<tag.GetPacketId()<<" received in RM time="<<Now().GetMilliSeconds()<<std::endl;
	NodeInfoHeader head;
	pkt->PeekHeader(head);
	Vector pos_d = head.getPosition();
	Ptr<Node> src = dev->GetNode();
	Ptr<MobilityModel> model = src->GetObject<MobilityModel> ();
	Vector src_position = model->GetPosition();
	if(CalculateDistance (pos_d, src_position) < Device_Transmission_Range){
		uint32_t src_id = src->GetId();
		if(isCluster[src_id]){
			if(head.getRiskFactor()<risk_factor[src_id]){
				isCluster[src_id] = false;
				std::map<uint32_t, std::vector<struct Cdq>*>::iterator j;
				j=CDQ.find(src_id);
				if(j!=CDQ.end()){
					j->second->clear();
					j->second->shrink_to_fit();
					delete j->second;
					j->second = 0;
					CDQ.erase(j);
				}
			}else{
				std::map<uint32_t, std::vector<struct Cdq>*>::iterator j;
				j=CDQ.find(src_id);
				Cdq cdq;
				cdq.send_id = send_id;
				cdq.comend=head.getRiskFactor();
				if(j==CDQ.end()){
					std::vector<Cdq> *vector = new std::vector<Cdq>;
					(*vector).push_back(cdq);
					cdq.send_id = src_id;
					cdq.comend=risk_factor[src_id];
					(*vector).push_back(cdq);
					CDQ.insert(std::make_pair(src_id,vector));
				}
				else{
					j->second->push_back(cdq);
				}
			}
		}
	}
}
void
MultipleChannelsExperiment::CdqPacketEvent(Ptr<NetDevice> dev, Ptr<const Packet> pkt, const Address &sender){
	StatsTag tag;
	bool result;
	result = pkt->FindFirstMatchingByteTag (tag);
	if (!result){
		NS_FATAL_ERROR ("the packet here shall have a stats tag");
	}
	  //std::cout<<"packet :"<<tag.GetPacketId()<<" received in CDQ time="<<Now().GetMilliSeconds()<<std::endl;

	NodeInfoHeader node_info_header;
	pkt->PeekHeader(node_info_header);
	CdqHeader header;
	pkt->PeekHeader(header);

    	Ptr<Node> src = dev->GetNode();
    	uint32_t src_id =src->GetId();
    	Ptr<MobilityModel> model = src->GetObject<MobilityModel> ();
    	Vector src_position = model->GetPosition();
    Vector pos_d = node_info_header.getPosition();
    bool flag = false;
	std::vector<struct Cdq> vec=header.getCdqVector();
	for(uint32_t i=0;i<vec.size();i++){
		Cdq cdq = vec.at(i);
		if(cdq.send_id==src_id)
			flag=true;
	}
	if(flag && (!isCluster[src_id])){
		std::map<uint32_t, std::vector<struct Cdq>*>::iterator j;
		j=CDQ1.find(src_id);
		//if the received invitations that from clusters more than one,choose the closest cluster to the origin
		if(j==CDQ1.end()){
			std::vector<struct Cdq>* vector = new std::vector<struct Cdq>;
			*vector = header.getCdqVector();
			CDQ1.insert(std::make_pair(src_id,vector));
			channel_table[src_id]=header.getChannel();
			nodes_cluster_address[src_id]=sender;
			nodes_cluster[src_id]=tag.getNodeId();
			isSuccefulCluster[src_id]=true;
			distant_from_cluster[src_id]=CalculateDistance(pos_d, src_position);
		}else{
			if(CalculateDistance (pos_d, src_position)<distant_from_cluster[src_id]){
				*(j->second)=header.getCdqVector();
				channel_table[src_id]=header.getChannel();
				nodes_cluster_address[src_id]=sender;
				nodes_cluster[src_id]=tag.getNodeId();
				isSuccefulCluster[src_id]=true;
				distant_from_cluster[src_id]=CalculateDistance(pos_d, src_position);
			}
		}
	}
}
void
MultipleChannelsExperiment::AmPacketEvent(Ptr<NetDevice> dev, Ptr<const Packet> pkt, const Address &sender){
	StatsTag tag;
	bool result;
	result = pkt->FindFirstMatchingByteTag (tag);
	if (!result){
		NS_FATAL_ERROR ("the packet here shall have a stats tag");
	}
	  //std::cout<<"packet :"<<tag.GetPacketId()<<" received in AM time="<<Now().GetMilliSeconds()<<std::endl;

	uint32_t src_id = dev->GetNode()->GetId();
	std::map<uint32_t, std::vector<uint32_t>*>::iterator j;
	j=CLUSTER.find(src_id);
	if(j==CLUSTER.end()){
		std::vector<uint32_t> * vector = new std::vector<uint32_t>;
		vector->push_back(tag.getNodeId());
		CLUSTER.insert(make_pair(src_id,vector));
	}else{
		j->second->push_back(tag.getNodeId());
	}
	//std::cout<<"I'm node "<<dev->GetNode()->GetId()<<" the "<<tag.getNodeId()<<" is my member!"<<std::endl;
}
void
MultipleChannelsExperiment::SchPacketEvent(Ptr<NetDevice> dev, Ptr<const Packet> pkt, const Address &sender){
	StatsTag tag;
	bool result;
	result = pkt->FindFirstMatchingByteTag (tag);
	if (!result){
		NS_FATAL_ERROR ("the packet here shall have a stats tag");
	}
	delaySum += (Now() - tag.GetSendTime()).GetMicroSeconds ();
	receiveSum++;
	std::map<uint32_t, std::vector<uint32_t>*>::iterator i;
	i=broadcastPackets.find(tag.getNodeId());
	if(i!=broadcastPackets.end()){
		std::vector<uint32_t>::iterator it;
		for(it=i->second->begin();it!=i->second->end();)
		{
			if(*it==dev->GetNode()->GetId())
				it=i->second->erase(it);
			else it++;
		}
		if(i->second->size()==1 && (i->second->at(0)==tag.getNodeId()))
			receives++;
	}

}


bool
MultipleChannelsExperiment::Receive (Ptr<NetDevice> dev, Ptr<const Packet> pkt, uint16_t mode, const Address &sender)
{
	  NS_LOG_FUNCTION (this << dev << pkt << mode << sender);
	  NS_ASSERT (mode == Packet_Number);
	  StatsTag tag;
	  bool result;
	  result = pkt->FindFirstMatchingByteTag (tag);
	  if (!result){
		  NS_FATAL_ERROR ("the packet here shall have a stats tag");
		  return false;
	  }
	  //std::cout<<"node "<<tag.GetPacketId()<<" Receive time"<<Seconds(getCurrentTime()-time_start).GetMilliSeconds()<<std::endl;
	  uint32_t packetType = tag.getPacketType();
	  switch (packetType){
	  case RmPacket:
		  RmPacketEvent(dev,pkt,tag.getNodeId());
		  break;
	  case CdqPacket:
		  CdqPacketEvent(dev,pkt,sender);
		  break;
	  case AmPacket:
		  AmPacketEvent(dev,pkt,sender);
		  break;
	  case SchPacket:
		  SchPacketEvent(dev,pkt,sender);
		  break;
	  default:
		  break;
	  }

	  return true;
}


void
MultipleChannelsExperiment::SendRmPackets (Time time, uint32_t channelNumber)
{
	NetDeviceContainer::Iterator i;
	bool flag=true;
	for (i = devices.Begin (); i != devices.End (); ++i)
	    {
	      Ptr<WaveNetDevice> sender = DynamicCast<WaveNetDevice> (*i);
	      double t = getCurrentTime()-time_start;
	      Simulator::Schedule(Seconds(rng->GetValue (t, t+0.020)),
	      				  	&MultipleChannelsExperiment::SendRm, this,
							sender, CCH,flag);
	  	flag=false;
	    }
}
void
MultipleChannelsExperiment::SendCdqPackets (uint32_t channelNumber)
{

	NetDeviceContainer::Iterator i;
		for (i = devices.Begin (); i != devices.End (); ++i)
		{
			Ptr<WaveNetDevice> sender = DynamicCast<WaveNetDevice> (*i);
			double t = getCurrentTime()-time_start;
			Simulator::Schedule(Seconds(rng->GetValue (t, t + 0.010)),
							&MultipleChannelsExperiment::SendCdq, this,
							sender);
		}

}
void
MultipleChannelsExperiment::SendAmPackets (Ptr<WaveNetDevice> sender, uint32_t channelNumber,bool falg)
{
	if(falg){
		for(uint32_t i=0;i<Car_number;i++){
			if(isCluster[i]){
				std::map<uint32_t, std::vector<struct Cdq>*>::iterator j;
				j=CDQ.find(i);
				if(j!=CDQ.end()){
					std::vector<struct Cdq>* vec = new std::vector<Cdq>;
					*vec=*(j->second);
					CDQ1.insert(std::make_pair(i,vec));
				}
			}
		}
	}
	  Ptr<Node> src = sender->GetNode();
		if(isSuccefulCluster[src->GetId()]){
			Time now = Now ();
			Ptr<Packet> packet = Create<Packet> (rm_size);
			StatsTag tag = StatsTag (++sequence, now,AmPacket,src->GetId());
			packet->AddByteTag (tag);
			Address dest = nodes_cluster_address[src->GetId()];
			bool result = false;
			TxInfo info = TxInfo (channelNumber);
			result = sender->SendX (packet, dest, Packet_Number, info);
	}
		//else
		//std::cout<<"channel_table[src->GetId()]= "<<channel_table[src->GetId()]<<std::endl;
}

void
MultipleChannelsExperiment::SendSchPackets (bool flag)
{
	if(flag){
		//std::cout<<Seconds(getCurrentTime()-time_start).GetMilliSeconds()<<" is in Sch interval!"<<std::endl<<std::endl<<std::endl;
		std::map<uint32_t, std::vector<struct Cdq>*>::iterator i;
		for (i = CDQ.begin(); i != CDQ.end();++i)
		{
			i->second->clear();
			i->second->shrink_to_fit();
			delete i->second;
			i->second = 0;
		}
		CDQ.clear();
		std::map<uint32_t, std::vector<uint32_t>*>::iterator ite;
		for (ite = CLUSTER.begin(); ite != CLUSTER.end();++ite)
		{
			std::cout<<"Cluster: "<<ite->first<<std::endl<<"The members: ";
			isSuccefulCluster[ite->first]=true;//簇头也是成功分簇的点，在这加上
			std::vector<uint32_t> *vector = ite->second;
			vector->push_back(ite->first);
			for(uint32_t i=0;i<vector->size();i++){
				std::cout<<vector->at(i)<<" ";
			}
			std::cout<<std::endl;
		}
		for (ite = broadcastPackets.begin(); ite != broadcastPackets.end();++ite)
		{
			ite->second->clear();
			ite->second->shrink_to_fit();
		}

		for(uint32_t i=0;i<Car_number;i++){
			if(isSuccefulCluster[i]){
				ite = CLUSTER.find(nodes_cluster[i]);
				if(ite!=CLUSTER.end()){
					std::vector<uint32_t> vec = *(ite->second);
					ite = broadcastPackets.find(i);
					if(ite!=broadcastPackets.end())
						*(ite->second) = vec;
				}
			}
		}
		for(int i=0;i<Car_number;i++){
			isCluster[i]=true;
		}
		for(int i=0;i<Car_number;i++){
			if(isSuccefulCluster[i])
				success_cluster++;
			else
				unsuccess_cluster++;
			isSuccefulCluster[i]=false;
		}
		for(int i=0;i<Car_number;i++){
			distant_from_cluster[i]=0.0;
		}
		for(int i=0;i<Car_number;i++){
		    	double risk=rng->GetValue (0.0, 10.0);
		    	risk_factor[i] = risk;
		 }
	}
	NetDeviceContainer::Iterator i;
	for (i = devices.Begin (); i != devices.End (); ++i){
		Ptr<WaveNetDevice> sender = DynamicCast<WaveNetDevice> (*i);
		uint32_t src_id = sender->GetNode()->GetId();
		std::map<uint32_t, std::vector<struct Cdq>*>::iterator i;
		i = CDQ1.find(src_id);
		if(i!=CDQ1.end()){
			//std::cout<<"node "<<src_id<<"Find in CDQ1"<<std::endl;
			std::vector<Cdq> *vector = i->second;
			uint32_t send_id=0;
			double comand=0.0;
			uint32_t index=0;
			for(uint32_t j=0;j<vector->size();j++){
				struct Cdq cdq=vector->at(j);
				if(cdq.comend>comand){
					send_id = cdq.send_id;
					comand = cdq.comend;
					index = j;
				}
			}
			if(comand>0){
				i->second->at(index).comend = 0.0;
				if(src_id == send_id){
					Time now = Now ();
					Ptr<Packet> packet = Create<Packet> (rm_size);
					StatsTag tag = StatsTag (++sequence, now,SchPacket,src_id);
					packet->AddByteTag (tag);
					Address dest = Mac48Address::GetBroadcast ();
					bool result = false;
					TxInfo info = TxInfo (channel_table[src_id]);
					result = sender->SendX (packet, dest, Packet_Number, info);
					if(result){
						sends++;
					}
					else
						std::cout<<"node "<<src_id<<"send faild in Sch! time= "<<Now().GetMilliSeconds()<<" IsSchInterval "<<sender->GetChannelCoordinator()->IsSchInterval()
						<<" IsGuardInterval "<<sender->GetChannelCoordinator()->IsGuardInterval()<<" IsCchInterval "<<sender->GetChannelCoordinator()->IsCchInterval()<<std::endl;		}
			}
		}
	}
}

void
MultipleChannelsExperiment::SendCdq(Ptr<WaveNetDevice> sender)
{
	uint32_t node_id=sender->GetNode()->GetId();
	if(isCluster[node_id]){
		uint32_t chanel = GetChannel();
		channel_table[node_id]=chanel;
		Time now = Now ();
		Ptr<Packet> packet = Create<Packet> (size);
		StatsTag tag = StatsTag (++sequence, now,CdqPacket,node_id);
		packet->AddByteTag (tag);
		Address dest = Mac48Address::GetBroadcast ();
		NodeInfoHeader node_info_header;
		CdqHeader header;
		std::map<uint32_t,std::vector<struct Cdq>*>::iterator j;
		j=CDQ.find(node_id);
		if(j!=CDQ.end()){
			std::vector<struct Cdq>* cdq_vector=j->second;
			//std::cout<<"Cluster: "<<node_id<<std::endl<<"the members: ";
			//for(int i=0;i<cdq_vector->size();i++){
			//	struct Cdq cdq=cdq_vector->at(i);
			//	std::cout<<cdq.send_id<<" ,";
			//}
			std::cout<<std::endl;
			Ptr<Node> src = sender->GetNode();
			Ptr<MobilityModel> model_src = src->GetObject<MobilityModel> ();
			const Vector pos_src = model_src->GetPosition ();
			node_info_header.setPosition(pos_src);
			node_info_header.setRiskFactor(risk_factor[node_id]);
			node_info_header.setSpeed(model_src->GetVelocity().x);
			packet->AddHeader(node_info_header);
			header.setChannel(chanel);
			header.setSize(cdq_vector->size());
			header.setCdqVector((*cdq_vector));
			packet->AddHeader(header);
			bool result = false;
			TxInfo info = TxInfo (CCH);
			result = sender->SendX (packet, dest, Packet_Number, info);
			if(!result){
				std::cout<<"CDQ发送失败"<<std::endl;
			}
		}else{
			isCluster[node_id]=false;
			std::cout<<"没找到CDQ"<<std::endl;
		}
	}
}
void MultipleChannelsExperiment::SendRm(Ptr<WaveNetDevice> sender,uint32_t channelNumber,bool flag)
{
	if(flag){
		for(int i=0;i< Car_number;i++){
			channel_table[i]=172;
		}
		std::map<uint32_t, std::vector<struct Cdq>*>::iterator ite;
		for (ite = CDQ1.begin(); ite != CDQ1.end();++ite)
		{
			ite->second->clear();
			delete ite->second;
			ite->second = 0;
		}
		CDQ1.clear();

		std::map<uint32_t, std::vector<uint32_t>*>::iterator ite_cluster;
		for (ite_cluster = CLUSTER.begin(); ite_cluster != CLUSTER.end();++ite_cluster)
		{
			ite_cluster->second->clear();
		    	ite_cluster->second->shrink_to_fit();
		  	delete ite_cluster->second;
		  	ite_cluster->second = 0;
		}
		CLUSTER.clear();
	}
		Time now = Now ();
		Ptr<Packet> packet = Create<Packet> (rm_size);
		StatsTag tag = StatsTag (++sequence, now,RmPacket,sender->GetNode()->GetId());
		packet->AddByteTag (tag);
		Address dest = Mac48Address::GetBroadcast ();
		Ptr<Node> src = sender->GetNode();
		Ptr<MobilityModel> model_src = src->GetObject<MobilityModel> ();
		const Vector pos_src = model_src->GetPosition ();
		NodeInfoHeader header;
		header.setPosition(pos_src);
		header.setRiskFactor(risk_factor[src->GetId()]);
		header.setSpeed(model_src->GetVelocity().x);
		packet->AddHeader(header);
    //std::cout<<"node "<<tag.GetPacketId()<<" send time"<<Seconds(getCurrentTime()-time_start).GetMilliSeconds()<<std::endl;
		bool result = false;
		TxInfo info = TxInfo (channelNumber);
		result = sender->SendX (packet, dest, Packet_Number, info);
		if(result){
			if(!(sender->GetChannelCoordinator()->IsRMInterval(now))){
				//std::cout<<"packet: "<<tag.GetPacketId()<<" is not in RM! time="<<now.GetMilliSeconds()<<std::endl;
			}
		}else
			std::cout<<"packet: "<<tag.GetPacketId()<<"send faild in RM!"<<std::endl;

}


uint32_t
MultipleChannelsExperiment::GetChannel()
{
	uint32_t channel = channel_current++;
	if (channel_current > 5)
		channel_current = channel_current % 6;
	return Channel[channel];
}

void
MultipleChannelsExperiment::InstallApplicationA (void)
{
  NS_LOG_FUNCTION (this);
  time_start = getCurrentTime();
  int flag = 0;
  int rm_flag = 0;
  int cdq_flag = 0;
  Ptr<WaveNetDevice> sender = DynamicCast<WaveNetDevice>(*(devices.Begin()));
  Ptr<ChannelCoordinator> coordinator_flag = sender->GetChannelCoordinator();
  NetDeviceContainer::Iterator i;
  while(getCurrentTime()-time_start < simulationTime)
  {
	  double now=getCurrentTime() - time_start;
	  Time t = Seconds(now);
	  if(coordinator_flag->IsCchInterval(t)&&flag==0){
		  if(coordinator_flag->IsRMInterval(t)&&rm_flag==0){
			  //std::cout<<t.GetMilliSeconds()<<" is in RM interval!"<<std::endl;
			  rm_flag=1;
			  SendRmPackets(t,CCH);
		  }
		  if(coordinator_flag->IsCDQInterval(t)&&cdq_flag==0){
			  cdq_flag = 1;
			  //std::cout<<t.GetMilliSeconds()<<" is in CDQ interval!"<<std::endl;
			 SendCdqPackets(CCH);
		  }
		  if(coordinator_flag->IsAMInterval(t)){
			  flag=1;
			  //std::cout<<t.GetMilliSeconds()<<" is in AM interval!"<<std::endl;
			  NetDeviceContainer::Iterator i;
			  bool flag_first = true;
			  for (i = devices.Begin (); i != devices.End (); ++i){
				  Ptr<WaveNetDevice> sender = DynamicCast<WaveNetDevice> (*i);
				  double t = getCurrentTime() - time_start;
			  	  Simulator::Schedule(Seconds(rng->GetValue (t, t + 0.008)),
			  			  		&MultipleChannelsExperiment::SendAmPackets, this,
			  					sender,CCH,flag_first);
			  	flag_first = false;
			  	}
		  }
	  }
	  else if(coordinator_flag->IsSchInterval(t) && flag==1){
		  //std::cout<<t.GetMilliSeconds()<<" is in SCH interval!"<<std::endl<<std::endl<<std::endl;
		  flag=0;
		  rm_flag=0;
		  cdq_flag=0;
		  Simulator::Schedule (Seconds (now-0.002), &MultipleChannelsExperiment::StartSch,this);
		  Simulator::Schedule (Seconds (now+0.05), &MultipleChannelsExperiment::StopSch,this);
		  bool temp = true;
		  double interval = coordinator_flag->GetSchInterval().GetSeconds()/Max_timeslot;
		  for(uint32_t N_timeslot=0;N_timeslot < Max_timeslot;N_timeslot++){
				  Simulator::Schedule(Seconds(rng->GetValue(now+interval*N_timeslot,now+interval*N_timeslot+(interval/3)*2)),
				  		  			  				  	  &MultipleChannelsExperiment::SendSchPackets, this,
				  										  temp);
			 temp=false;
		  }
	  }

  }
}



void
MultipleChannelsExperiment::Run (void)
{
  NS_LOG_FUNCTION (this);
  NS_LOG_DEBUG ("simulation configuration arguments: ");
  NS_LOG_UNCOND ("      PDR  AverageDelay SystemThroughput  AverageThroughput");
  {
    NS_LOG_DEBUG ("configuration A:");
    {
    pdr = delay = system_throughput = average_throughput = 0.0;
    for (uint32_t r = 0; r != run; r++)
    {
    	RngSeedManager::SetSeed (1);
    	RngSeedManager::SetRun (17+r);
    	InitStats ();
    	CreateNodes ();
    	SetupMobility ();
    	CreateWaveDevice ();
    	InstallApplicationA ();
    	Simulator::Stop (Seconds (simulationTime+1));
    	AnimationInterface anim ("wave-Culstering.xml"); // Mandatory
    	//anim.SetMobilityPollInterval (Seconds (1));
    	Simulator::Run ();
    	Simulator::Destroy ();
    	Stats (17+r);
    	Destroy();
    }
    NS_LOG_UNCOND ("A :  " << pdr/run<< " " << delay/run << " " << system_throughput/run << " " << average_throughput/run);
    }
  }
}


void
MultipleChannelsExperiment::InitStats (void)
{
  // used for sending packets randomly
  rng = CreateObject<UniformRandomVariable> ();
  rng->SetStream (1);
  for(uint32_t i=0;i<Car_number;i++){
	  std::vector<uint32_t> *vec = new std::vector<uint32_t>;
	  broadcastPackets.insert(std::make_pair(i,vec));
  }

  std::map<uint32_t, std::vector<struct Cdq>*>::iterator ite;
  			for (ite = CDQ.begin(); ite != CDQ.end();++ite)
  			{
  				ite->second->clear();
  				ite->second->shrink_to_fit();
  				delete ite->second;
  				ite->second = 0;
  			}
  			CDQ.clear();

  std::map<uint32_t, std::vector<struct Cdq>*>::iterator ite_cdq1;
  for (ite_cdq1 = CDQ1.begin(); ite_cdq1 != CDQ1.end();++ite_cdq1)
  {
	  ite_cdq1->second->clear();
	  ite_cdq1->second->shrink_to_fit();
	  delete ite_cdq1->second;
	  ite_cdq1->second = 0;
  }
  CDQ1.clear();
  std::map<uint32_t, std::vector<uint32_t>*>::iterator ite_cluster;
    for (ite_cluster = CLUSTER.begin(); ite_cluster != CLUSTER.end();++ite_cluster)
    {
    	ite_cluster->second->clear();
    	ite_cluster->second->shrink_to_fit();
  	delete ite_cluster->second;
  	ite_cluster->second = 0;
    }
    CLUSTER.clear();

    for(int i=0;i<Car_number;i++){
    		isCluster[i]=true;
    	}
    for(int i=0;i<Car_number;i++){
    		isSuccefulCluster[i]=false;
    	}
    for(int i=0;i<Car_number;i++){
    		distant_from_cluster[i]=true;
   	}
    for(int i=0;i<Car_number;i++){
    		risk_factor[i]=true;
    }
    for(int i=0;i<Car_number;i++){
    		double risk=rng->GetValue (0.0, 10.0);
    	    	risk_factor[i] = risk;
    }

  receivedPackets.clear();
  success_cluster=0;
  unsuccess_cluster=0;
  sends = 0;
  sequence = 0;
  receives = 0;
  delaySum = 0;
  receiveSum = 0;

  //Simulator::ScheduleDestroy (&MultipleChannelsExperiment::StatQueuedPackets, this);
}

void
MultipleChannelsExperiment::Destroy (void){
	std::map<uint32_t, std::vector<uint32_t> *>::iterator i;
	for (i = broadcastPackets.begin(); i != broadcastPackets.end();++i)
	{
		i->second->clear();
		i->second->shrink_to_fit();
		delete i->second;
		i->second = 0;
	}
	broadcastPackets.clear();

	  std::map<uint32_t, std::vector<uint32_t>*>::iterator ite_cluster;
	    for (ite_cluster = CLUSTER.begin(); ite_cluster != CLUSTER.end();++ite_cluster)
	    {
	    	ite_cluster->second->clear();
	    	ite_cluster->second->shrink_to_fit();
	  	delete ite_cluster->second;
	  	ite_cluster->second = 0;
	    }
	    CLUSTER.clear();

	std::map<uint32_t, std::vector<struct Cdq>*>::iterator ite;
			for (ite = CDQ.begin(); ite != CDQ.end();++ite)
			{
				ite->second->clear();
				ite->second->shrink_to_fit();
				delete ite->second;
				ite->second = 0;
			}
			CDQ.clear();

			std::map<uint32_t, std::vector<struct Cdq>*>::iterator ite_cdq1;
			  for (ite_cdq1 = CDQ1.begin(); ite_cdq1 != CDQ1.end();++ite_cdq1)
			  {
				  ite_cdq1->second->clear();
				  ite_cdq1->second->shrink_to_fit();
				  delete ite_cdq1->second;
				  ite_cdq1->second = 0;
			  }
			  CDQ1.clear();

}

void
MultipleChannelsExperiment::Stats (uint32_t randomNumber)
{
  NS_LOG_FUNCTION (this);
  NS_LOG_DEBUG (" simulation result:(run number = " << randomNumber << ")");
  NS_LOG_DEBUG (" nodes = " << nodesNum);
  NS_LOG_DEBUG (" simulation time = " << simulationTime << "s");
  NS_LOG_DEBUG (" AverageClustering = " << (double)success_cluster/(double)(success_cluster+unsuccess_cluster)*100 << "%");
  NS_LOG_DEBUG (" send packets = " << sends);
  NS_LOG_DEBUG (" receives packets = " << receives);
  NS_LOG_DEBUG (" receiveSum packets = " << receiveSum);
  NS_LOG_DEBUG (" lost packets = " << (sends - receives));
  NS_LOG_DEBUG (" the sum delay of all received packets = " << delaySum << "micros");


  // second show performance result
  // stats PDR (packet delivery ratio)
  double PDR = receives / (double)(sends);
  // stats average delay
  double AverageDelay = delaySum / receiveSum / 1000.0;
  // stats system throughput, Mbps
  double Throughput = receives * size * 8 / simulationTime / 1000.0 / 1000.0;
  // stats average throughput, kbps
  double AverageThroughput = Throughput / (double)nodesNum * 1000;

  NS_LOG_DEBUG (" PDR = " << PDR);
  NS_LOG_DEBUG (" AverageDelay = " << AverageDelay << "ms");
  NS_LOG_DEBUG (" Throughput = " << Throughput << "Mbps");
  NS_LOG_DEBUG (" AverageThroughput = " << AverageThroughput << "Kbps");

  pdr+=PDR;
  delay+=AverageDelay;
  system_throughput+=Throughput;
  average_throughput+=AverageThroughput;

}
void
MultipleChannelsExperiment::StartSch ()
{
	NetDeviceContainer::Iterator it;
	for(it = devices.Begin ();it!=devices.End();it++){
		Ptr<WaveNetDevice> sender = DynamicCast<WaveNetDevice> (*it);
		const SchInfo schInfo = SchInfo (channel_table[sender->GetNode()->GetId()], false, EXTENDED_ALTERNATING);
		sender->StartSch(schInfo);
	}

}
void
MultipleChannelsExperiment::StopSch ()
{
	NetDeviceContainer::Iterator it;
	for(it = devices.Begin ();it!=devices.End();it++){
		Ptr<WaveNetDevice> sender = DynamicCast<WaveNetDevice> (*it);
		const SchInfo schInfo = SchInfo (channel_table[sender->GetNode()->GetId()], false, EXTENDED_ALTERNATING);
		sender->StopSch(channel_table[sender->GetNode()->GetId()]);
	}

}

int
main (int argc, char *argv[])
{
  LogComponentEnable ("WaveMultipleChannel", LOG_DEBUG);

  MultipleChannelsExperiment experiment;
  if (experiment.Configure (argc, argv))
    {
      experiment.Run ();
    }
  else
    {
      experiment.Usage ();
    }

  return 0;
}
