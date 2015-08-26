package com.epcc.monc;

import java.util.HashMap;
import java.util.Map;

public class CompilerTests {
	private String name;
	private Map<String, Boolean> testRuns;
	private Map<String, StringBuilder> testRunMessages;
	private Map<String, Integer> testCases;
	private Map<String, Integer> testCasePasses;
	private boolean pass;
	private boolean passCompile;
	private StringBuilder compileMessages;
	
	public CompilerTests(String name) {
		this.name = name;
		this.testRuns = new HashMap<>();
		this.testRunMessages = new HashMap<>();
		this.testCases = new HashMap<>();
		this.testCasePasses = new HashMap<>();
		this.pass = true;
		this.passCompile = true;
		compileMessages = new StringBuilder();
	}
	
	public void setPass(boolean pass) {
		this.pass = pass;
	}
	
	public void setPassCompile(boolean passCompile) {
		this.passCompile = passCompile;
	}
	
	public boolean getPassCompile() {
		return this.passCompile;
	}
	
	public boolean getPass() {
		return this.pass;
	}
	
	public void addCompileMessage(String message) {
		compileMessages.append(message).append("<br>");
	}
	
	public String getCompileMessage() {
		return compileMessages.toString();
	}
	
	public void setTestCasePasses(String name, int passes) {
		testCasePasses.put(name, passes);
	}
	
	public void setTestCases(String name, int cases) {
		testCases.put(name, cases);
	}
	
	public int getTotalTestCases() {
		int cases = 0;
		for (int singleCases : testCases.values()) {
			cases+=singleCases;
		}
		return cases+1;
	}
	
	public int getTotalTestCasePasses() {
		int passes = passCompile ? 1 : 0;
		for (int singleCases : testCasePasses.values()) {
			passes+=singleCases;
		}
		return passes;
	}
	
	public int getTestCasePasses(String name) {
		if (testCasePasses.containsKey(name)) {
			return testCasePasses.get(name); 
		}
		return 0;
	}
	
	public int getTestCases(String name) {
		if (testCases.containsKey(name)) {
			return testCases.get(name); 
		}
		return 0;
	}
	
	public void addTest(String name, boolean pass) {
		testRuns.put(name, pass);
	}
	
	public void addTestMessage(String name, String messageLine) {
		if (!testRunMessages.containsKey(name)) testRunMessages.put(name, new StringBuilder());
		testRunMessages.get(name).append(messageLine).append("<br>");
	}
	
	public String getTestMessage(String name) {
		if (testRunMessages.containsKey(name)) {
			return testRunMessages.get(name).toString();
		} else {
			return null;
		}
	}
	
	public Map<String, Boolean> getTestRuns() {
		return testRuns;
	}
	
	public boolean getTestSuccess(String name) {
		return testRuns.get(name);
	}
	
	public String getName() {
		return name;
	}
}
