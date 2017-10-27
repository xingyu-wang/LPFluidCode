#include "boundary.h"
#include <iostream>


////////////////////////////////////////////////////////////////////////////////////////
// Start of BoundaryFactory
////////////////////////////////////////////////////////////////////////////////////////

BoundaryFactory& BoundaryFactory::instance() { 
	static BoundaryFactory bFactory;
	return bFactory;
}

void BoundaryFactory::registerBoundary(std::string name, GeoCreateFunc func) {
	
	bTable.insert({name,func});

}

Boundary* BoundaryFactory::createBoundary(std::string name) {
	const auto result = bTable.find(name);
	if(result==bTable.end()) { // the goemetry class name is not registered
		std::cout<<"This boundary class name is not registered!!!"<<std::endl;
		return nullptr;
	}
	return (result->second)();
}

////////////////////////////////////////////////////////////////////////////////////////
// End of BoundaryFactory
////////////////////////////////////////////////////////////////////////////////////////

