package org.mg.cdklib.match;

import java.util.HashMap;
import java.util.Map;

import org.mg.javalib.util.StringUtil;
import org.openscience.cdk.fingerprint.SubstructureFingerprinter;

public class SubstructureCountMatcherProvider
{
	public static final String FOLDER = System.getProperty("user.home") + "/results/cdk-match/";

	public static SubstructureCountMatcher getInstance(SubstructureFingerprinter sub)
	{
		String key = sub.getClass().getSimpleName();
		if (instances.containsKey(key))
			return instances.get(key);

		String smarts[] = new String[sub.getSize()];
		for (int i = 0; i < sub.getSize(); i++)
			smarts[i] = sub.getSubstructure(i);
		return getInstance(smarts, key);
	}

	public static SubstructureCountMatcher getInstance(String smarts[])
	{
		StringBuffer b = new StringBuffer();
		for (String smrt : smarts)
		{
			b.append(smrt);
			b.append('\n');
		}
		return getInstance(smarts, StringUtil.getMD5(b.toString()));
	}

	private static SubstructureCountMatcher getInstance(String smarts[], String key)
	{
		if (!instances.containsKey(key))
		{
			String file = FOLDER + DefaultSubstructureCountMatcher.getVersion() + "#" + key;
			instances.put(key, new PersistentSubstructureCountMatcher(smarts, file));
		}
		return instances.get(key);
	}

	static Map<String, SubstructureCountMatcher> instances = new HashMap<>();

}