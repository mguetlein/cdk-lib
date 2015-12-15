package org.mg.cdklib.cfp;

import java.io.Serializable;

/**
 * this class wraps fragments-id to avoid confusion with compound-idx in cfpminer
 * 
 * when folding is enabled, the id refers to the bitVector position, otherwise it is the CFP fragment id
 */
public class CFPFragment implements Serializable
{
	private static final long serialVersionUID = 1L;

	final private int id;

	public CFPFragment(int id)
	{
		this.id = id;
	}

	public int getId()
	{
		return id;
	}

	@Override
	public int hashCode()
	{
		return id;
	}

	@Override
	public boolean equals(Object obj)
	{
		return obj instanceof CFPFragment && id == ((CFPFragment) obj).id;
	}

	@Override
	public String toString()
	{
		return Integer.toString(id);
	}
}
