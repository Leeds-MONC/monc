package com.epcc.monc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class UnitParser {
	
	private static Map<String, CompilerTests> compilerTests;
	private static Map<String, String> unitTestURLs;
	private static Set<String> unitTestCategories;
	private static String runTimestamp;	
	private static String svnVersion;
	private static boolean overallPass;
	private static int totalCases, totalPasses;	
	
	static {
		compilerTests = new HashMap<>();
		unitTestURLs = new HashMap<>();
		unitTestCategories = new HashSet<>();
		overallPass = true;
	}
	
	public static void main(String args[]) {
		try {
			parse(args[0]);
			String html = HtmlProcessor.generate(runTimestamp, compilerTests, unitTestCategories, unitTestURLs, svnVersion, overallPass,
					totalCases, totalPasses);
			writeHTML(args[1], html);
		} catch (IOException e) {			
			e.printStackTrace();
		}
	}
	
	private static void writeHTML(String filename, String html) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
		bw.write(html);
		bw.close();
	}
	
	private static void parse(String filename) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(filename));
		String line, currentCompiler=null, currentTest=null;
		boolean messageMode=false, compileMode=false;
		while ((line = br.readLine()) != null) {
			if (line.startsWith("[Result]")) {
				String consituientParts[] = line.split(":");
				currentCompiler = consituientParts[0].substring(8);
				if (!compilerTests.containsKey(currentCompiler)) compilerTests.put(currentCompiler, new CompilerTests(currentCompiler));
				compilerTests.get(currentCompiler).addTest(consituientParts[1], consituientParts[3].equals("1"));
				unitTestCategories.add(consituientParts[1]);
				currentTest = consituientParts[1];
				if (!unitTestURLs.containsKey(currentTest)) unitTestURLs.put(currentTest, consituientParts[2]);
				messageMode=true;
			} else if (line.startsWith("[Compiler]")) {
				String consituientParts[] = line.split(":");
				currentCompiler = consituientParts[0].substring(10);
				if (!compilerTests.containsKey(currentCompiler)) compilerTests.put(currentCompiler, new CompilerTests(currentCompiler));
				compilerTests.get(currentCompiler).setPass(consituientParts[1].equals("Pass"));
			} else if (line.startsWith("[Entire Build]")) {
				String consituientParts[] = line.split(":");
				currentCompiler = consituientParts[0].substring(14);
				if (!compilerTests.containsKey(currentCompiler)) compilerTests.put(currentCompiler, new CompilerTests(currentCompiler));
				compilerTests.get(currentCompiler).setPassCompile(consituientParts[1].equals("Pass"));
				totalCases++;
				if (compilerTests.get(currentCompiler).getPassCompile()) totalPasses++;
				compileMode=true;
			} else if (line.startsWith("[SVN]")) {
				String consituientParts[] = line.split(":");
				svnVersion = consituientParts[1];
			} else if (line.startsWith("[Entire]")) {
				overallPass = line.substring(8).equals("Pass");
			} else if (line.startsWith("[Start]")) {
				runTimestamp = line.substring(7).replace('-', ' ');
			} else if (line.equals("--End Message--")) {
				messageMode = compileMode = false;
			} else if ((messageMode || compileMode) && currentCompiler != null && currentTest != null){
				if (messageMode) compilerTests.get(currentCompiler).addTestMessage(currentTest, line);
				if (compileMode) compilerTests.get(currentCompiler).addCompileMessage(line);
				if (line.trim().startsWith("Successful cases   / total cases")) {
					String consituents[] = line.split(":");
					String passAndFail[]=consituents[1].split("/");
					passAndFail[0] = passAndFail[0].replace('[', ' ').trim();
					passAndFail[1] = passAndFail[1].replace(']', ' ').trim();
					compilerTests.get(currentCompiler).setTestCasePasses(currentTest, Integer.valueOf(passAndFail[0]));
					compilerTests.get(currentCompiler).setTestCases(currentTest, Integer.valueOf(passAndFail[1]));
					totalCases+=Integer.valueOf(passAndFail[1]);
					totalPasses+=Integer.valueOf(passAndFail[0]);
				}
			}
		}
		br.close();
	}
}
