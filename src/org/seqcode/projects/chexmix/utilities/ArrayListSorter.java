package org.seqcode.projects.chexmix.utilities;

public class ArrayListSorter implements Comparable<ArrayListSorter>{
	
	public static boolean ASC = true;
	public static boolean DECS = false;
	private Object key;
	private float val;
	
	public ArrayListSorter(){}
	
	public ArrayListSorter(Object key , float val){
		super();
		this.key = key;
		this.val = val;
	}
	
	public Object getKey(){
		return key;
	}
	
	public float getVal(){ 
		return val;
	}

	public void setKey(Object key){
		this.key = key;
	}
	
	public void setVal(float val){
		this.val = val;
	}
	
	public int compareTo (ArrayListSorter compareObj){
		
		float compareVal = ((ArrayListSorter) compareObj).getVal();
		
		// descending order
		if (compareVal > this.val){
			return 1;
		}else{
			return -1;
		}
	}
}
