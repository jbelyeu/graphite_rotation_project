#ifndef GRAPHITE_IPATH_H
#define GRAPHITE_IPATH_H

#include "core/alignment/IAlignment.h"
#include "core/allele/IAllele.h"
#include "core/util/Noncopyable.hpp"


namespace graphite
{
	class IPath : Noncopyable
	{
	public:
		typedef std::shared_ptr< IPath > SharedPtr;
		IPath() {}
		virtual ~IPath() {}

		virtual std::vector< IAllele::SharedPtr > getAllelePath() = 0;
		virtual uint32_t getPathSWPercent() = 0;
		virtual IAlignment::SharedPtr getAlignment() = 0;
		virtual size_t getHash() = 0;

		virtual void addAlleleToPath(IAllele::SharedPtr allelePtr) = 0;
		virtual void setPathSWPercent(uint32_t swPercent) = 0;
		virtual void setAlignment(IAlignment::SharedPtr allelePtr) = 0;
	};
}

#endif //GRAPHITE_IPATH_H
