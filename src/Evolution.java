import java.util.concurrent.ThreadLocalRandom;

public class Evolution {
	
	private final double min = -100;
	private final double max = 100;
	private final double cr = 0.9;
	private final double f = 0.8;
	private final int d = 10;
	private final int np = 100;
	private final int termination = 3000 * d;
	private double[][] vector = new double[np][d];
	private double[] target = new double[d];
	private double[] mutant = new double[d];
	private double[] child = new double[d];
	private double bestFit;
	private int nfc = 0;
	private double optimal = 0;

	public void initializePopulation(){
		for(int i=0; i < np; i++){
			for(int j=0; j < d; j++){
				vector[i][j] = ThreadLocalRandom.current().nextDouble(min, max + 1);
			}
		}
		
//		for(int i=0; i<np; i++){
//			System.out.println(vector[i]);
//			for(int j=0; j<d; j++){
//				System.out.println(vector[i][j]);
//			}
//		}
		setTargetVector();
	}
	
	private void setTargetVector(){
		if(nfc <= termination){
		target = vector[0];
//		for(double num : target){
//			System.out.println("Target: " + num);
//		}
		beginMutation();
		} else {
			for(double num : target){
				System.out.println("Target: " + num);
			}
			
			System.out.println("Finished." + optimal);
		}
	}
	
	private void beginMutation(){
		double[] rand1 = vector[ThreadLocalRandom.current().nextInt(1, d)];
		double[] rand2 = vector[ThreadLocalRandom.current().nextInt(1, d)];
		double[] rand3 = vector[ThreadLocalRandom.current().nextInt(1, d)];
				
		for(int i=0; i<d; i++){
			mutant[i] = rand1[i] + f * (rand2[i] - rand3[i]);
		}
//		for(double m : mutant){
//			System.out.println("Mutant: " + m);
//		}
		determineCrossover();
	}
	
	private void determineCrossover(){
		for(int i=0; i<d; i++){
			double rand = ThreadLocalRandom.current().nextDouble(0, 1);
//			System.out.println(rand);
			if(rand < cr){
				child[i] = mutant[i];
			} else {
				child[i] = target[i];
			}
		}
//		for(double c : child){
//			System.out.println("Child: " + c);
//		}
		testFunc1();
	}
	
	private void testFunc1(){
		//High Conditioned Elliptic Function
		double targetTot = 0;
		for(int i=0; i<d; i++){
			targetTot += Math.pow(Math.pow(10, 6), ((i-1)/d-1)) * Math.pow(target[i], 2);
		}
		
		double childTot = 0;
		for(int i=0; i<d; i++){
			childTot += Math.pow(Math.pow(10, 6), ((i-1)/d-1)) * Math.pow(child[i], 2);
		}
		
		check(targetTot, childTot);
	}
	
	private void testFunc2(){
		//Bent Cigar Function
		double targetTot = 0;
		for(int i=1; i<d; i++){
			targetTot += Math.pow(target[i], 2);
		}
		targetTot *= Math.pow(10, 6);
		targetTot = Math.pow(target[0], 2);
		
		double childTot = 0;
		for(int i=1; i<d; i++){
			childTot += Math.pow(child[i], 2);
		}
		childTot *= Math.pow(10, 6);
		childTot = Math.pow(child[0], 2);
		
		check(targetTot, childTot);
	}
	
	public void testFunc3(){
		//Discus Function
		double targetTot = 0;
		targetTot = Math.pow(10, 6) * Math.pow(target[0], 2);
		for(int i=1; i<d; i++){
			targetTot += Math.pow(target[i], 2);
		}
		
		double childTot = 0;
		childTot = Math.pow(10, 6) * Math.pow(child[0], 2);
		for(int i=1; i<d; i++){
			childTot += Math.pow(child[i], 2);
		}	
		check(targetTot, childTot);
	}
	
	public void testFunc4(){
		//Rosenbrock Function
		double targetTot = 0;
		for(int i=0; i<d-1; i++){
			targetTot += (100 * Math.pow(((Math.pow(target[i], 2)) - target[i+1]), 2) + Math.pow(target[i]-1, 2));
		}
		
		double childTot = 0;
		for(int i=0; i<d-1; i++){
			childTot += (100 * Math.pow(((Math.pow(child[i], 2)) - child[i+1]), 2) + Math.pow(child[i]-1, 2));
		}	
		check(targetTot, childTot);
	}
	
	public void testFunc5(){
		//Ackley's Function
		double targetTot1 = 0;
		for(int i=0; i<d; i++){
			targetTot1 += Math.pow(target[i], 2);
		}
		targetTot1 = targetTot1/d;
		targetTot1 *= -2;
		targetTot1 = Math.exp(targetTot1);
		targetTot1 *= 20;
		
		double targetTot2 = 0;
		for(int i=0; i<d; i++){
			targetTot2 += Math.cos(2*Math.PI*target[i]);
		}
		targetTot2 = targetTot2 / d;
		targetTot2 = Math.exp(targetTot2);
		
		double targetTot = targetTot1 - targetTot2 + 20 + Math.E;
		
		double childTot1 = 0;
		for(int i=0; i<d; i++){
			childTot1 += Math.pow(child[i], 2);
		}
		childTot1 = childTot1/d;
		childTot1 *= -2;
		childTot1 = Math.exp(childTot1);
		childTot1 *= 20;
		
		double childTot2 = 0;
		for(int i=0; i<d; i++){
			childTot2 += Math.cos(2*Math.PI*child[i]);
		}
		childTot2 = childTot2 / d;
		childTot2 = Math.exp(childTot2);
		
		double childTot = targetTot1 - targetTot2 + 20 + Math.E;
		
		check(targetTot, childTot);
	}
	
	private void testFunc6(){
		//Weierstrass Function
		double targetTot = 0;
		double targetTot1 = 0, targetTot2 = 0;
		for(int i=0; i<d; i++){
			for(int k=0; k<=20; k++){
				targetTot1 += Math.pow(0.5, k) * Math.cos(2*Math.PI*Math.pow(3, k)*(target[i]+0.5));
			}
			targetTot += targetTot1;
		}
		for(int k=0; k<=20; k++){
			targetTot2 += Math.pow(0.5, k) * Math.cos(2*Math.PI*Math.pow(3, k)*0.5);
		}
		targetTot2 *= d;
		targetTot -= targetTot2;
		
		double childTot = 0;
		double childTot1 = 0, childTot2 = 0;
		for(int i=0; i<d; i++){
			for(int k=0; k<=20; k++){
				childTot1 += Math.pow(0.5, k) * Math.cos(2*Math.PI*Math.pow(3, k)*(child[i]+0.5));
			}
			childTot += childTot1;
		}
		for(int k=0; k<=20; k++){
			childTot2 += Math.pow(0.5, k) * Math.cos(2*Math.PI*Math.pow(3, k)*0.5);
		}
		childTot2 *= d;
		childTot -= childTot2;
		
		check(targetTot, childTot);
	}
	
	public void testFunc7(){
		//Griewank Function
		double targetSum = 0;
		double targetProd = 1;
		for(int i=1; i<d+1; i++){
			targetSum += Math.pow(target[i-1], 2) / 4000;
			targetProd *= Math.cos(target[i-1]/Math.sqrt(i));
		}
		double targetTot = targetSum - targetProd + 1;
		
		double childSum = 0;
		double childProd = 1;
		for(int i=1; i<d+1; i++){
			childSum += Math.pow(child[i-1], 2) / 4000;
			childProd *= Math.cos(child[i-1]/Math.sqrt(i));
		}
		double childTot = childSum - childProd + 1;
		
		check(targetTot, childTot);
	}
	
	public void testFunc8(){
		//Rastrigin Function
		double targetTot = 0;
		for(int i=0; i<d; i++){
			targetTot += Math.pow(target[i], 2) - (10*Math.cos(2*Math.PI*target[i])) + 10;
		}
		
		double childTot = 0;
		for(int i=0; i<d; i++){
			childTot += Math.pow(child[i], 2) - (10*Math.cos(2*Math.PI*child[i])) + 10;
		}
		
		check(targetTot, childTot);
	}
	
	public void testFunc9(){
		//Katsuura Function
		double targetSum = 0;
		double targetProd = 1;
		for(int i=0; i<d; i++){
			for(int j=1; j<=32; j++){
				targetSum += Math.pow(Math.abs(Math.pow(2, j)*target[i] - Math.floor(Math.pow(2, j)*target[i])/Math.pow(2, j)), 10/Math.pow(d, 1.2));
			}
			targetProd *= targetSum * i + 1;
		}
		double targetTot = (10/Math.pow(d, 2)) * targetProd - (10/Math.pow(d, 2));
		
		double childSum = 0;
		double childProd = 1;
		for(int i=0; i<d; i++){
			for(int j=1; j<=32; j++){
				childSum += Math.pow(Math.abs(Math.pow(2, j)*child[i] - Math.floor(Math.pow(2, j)*child[i])/Math.pow(2, j)), 10/Math.pow(d, 1.2));
			}
			childProd *= childSum * i + 1;
		}
		double childTot = (10/Math.pow(d, 2)) * childProd - (10/Math.pow(d, 2));
		
		check(targetTot, childTot);
	}
	
	public void check(double t, double c){
		if(t > c){
			optimal = c;
			//System.out.println("--C Optimal" + optimal);
			vector[0] = child;
			nfc++;
			System.out.println("----NFC----:" + nfc);
			setTargetVector();
		} else {
			optimal = t;
			//System.out.println("--T Optimal" + optimal);
			vector[0] = target;
			nfc++;
			System.out.println("----NFC----:" + nfc);
			setTargetVector();
		}
	}
}

