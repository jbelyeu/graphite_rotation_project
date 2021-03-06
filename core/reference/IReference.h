#ifndef GRAPHITE_IREFERENCE_H
#define GRAPHITE_IREFERENCE_H

#include <stdint.h>
#include <memory>

#include "core/region/Region.h"
#include "core/util/Types.h"

namespace graphite
{

	class IReference : private Noncopyable
	{
	public:
		typedef std::shared_ptr<IReference> SharedPtr;

	    IReference(){}
		virtual ~IReference() {}

		virtual const char* getSequence() = 0;
		virtual size_t getSequenceSize() = 0;

		Region::SharedPtr getRegion() { return this->m_region; }

	protected:
		Region::SharedPtr m_region;

	};
}

#endif // GRAPHITE_IREFERENCE_H
