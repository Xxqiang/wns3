/*
 * test.cc
 *
 *  Created on: 2018年4月7日
 *      Author: qiangge666
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


int main(){
	 std::vector<int> vec;
	 int i;
	 i=1;
	 vec.push_back(i);
	 i=2;
	 vec.push_back(i);
	 std::cout<<"size= "<<vec.size()<<std::endl;
	 for(int j=0; j<vec.size();j++){
		std::cout<<vec.at(j)<<std::endl;
	 }
	return 0;
}

