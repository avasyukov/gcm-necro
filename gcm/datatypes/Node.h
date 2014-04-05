#ifndef _GCM_NODE_H
#define _GCM_NODE_H  1

#define PLACEMENT_TYPE_MASK 2
enum PlacementType
{
	Local = 2,
	Remote = 0
};

#define IS_USED_MASK 4

#define IS_BORDER_MASK 8

#define CONTACT_TYPE_MASK 1
enum ContactType
{
	Free = 0,
	InContact = 1
};

#define ENGINE_OWNERSHIP_MASK 48
enum EngineOwner
{
	GCM = 16,
	SPH = 32
};

#include "Basis.h"
#include "ContactData.h"
#include <assert.h>

//TODO: remove asserts from setXxx methods after some more tests
class Node
{
	friend class VTKSnapshotWriter;
	friend class TetrMesh_1stOrder;
public:
	int local_zone_num;
	int remote_zone_num;
	int local_num;
	int remote_num;
	int absolute_num;
	contact_state* contact_data;
	basis* local_basis;
	float coords[3];
	float fixed_coords[3];

	Node ()
	{
		node_flags = 0;
	}

	bool inline isInContact ()
	{
		return InContact == (node_flags & CONTACT_TYPE_MASK);
	}

	void inline setContactType (ContactType type)
	{
		node_flags &= (~CONTACT_TYPE_MASK);
		node_flags |= (CONTACT_TYPE_MASK & type);

		assert (InContact == type ? isInContact () : !isInContact ());
	}

	/**
	 *
	 * @return <code>true</code> if this Node is used and its placement is Local
	 */
	bool inline isLocal ()
	{
		return isUsed () && Local == (node_flags & PLACEMENT_TYPE_MASK);
	}

	/**
	 *
	 * @return <code>true</code> if this Node is used and its placement is Remote
	 */
	bool inline isRemote ()
	{
		return isUsed () && Remote == (node_flags & PLACEMENT_TYPE_MASK);
	}

	/**
	 * Set node placement to specified value and also marks node as used
	 * @param placement new node placement
	 */
	void inline setPlacement (PlacementType placement)
	{
		setUsed (true);

		node_flags &= (~PLACEMENT_TYPE_MASK);
		node_flags |= (PLACEMENT_TYPE_MASK & placement);

		assert (Local == placement ? isLocal () : isRemote ());
	}

	bool inline isUsed ()
	{
		return 0 != (node_flags & IS_USED_MASK);
	}

	void inline setUsed (bool used)
	{
		if (used) node_flags |= IS_USED_MASK;
		else node_flags &= (~IS_USED_MASK);

		assert (used ? isUsed () : !isUsed ());
	}

	void inline setIsBorder (bool border)
	{
		if (border) node_flags |= IS_BORDER_MASK;
		else node_flags &= (~IS_BORDER_MASK);

		assert (border ? isBorder () : !isBorder ());
	}

	bool inline isBorder ()
	{
		return 0 != (node_flags & IS_BORDER_MASK);
	}

	bool inline isInner ()
	{
		return !isBorder ();
	}

	bool inline isOwnedBy (EngineOwner owner)
	{
		assert (0 != (node_flags & ENGINE_OWNERSHIP_MASK));//Node should be owned at least by someone
		return 0 != (node_flags & ENGINE_OWNERSHIP_MASK & owner);
	}

	bool inline addOwner (EngineOwner owner)
	{
		node_flags |= (ENGINE_OWNERSHIP_MASK & owner);
	}
protected:
private:
	/**
	 * This method is only to be used to dump mesh state into a file or any other output stream
     * @return
     */
	unsigned int inline getFlags ()
	{
		return node_flags;
	}

	/**
	 * This method is only supposed to be used to read mesh state from a file or any other input stream
     * @param flags
     */
	void inline setFlags (unsigned int flags)
	{
		node_flags = flags;
	}

	unsigned int node_flags;
};

#endif
