import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import cern.jet.random.sampling.RandomSampler;
import cern.jet.random.Exponential;
import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;
import cern.jet.random.engine.RandomEngine;

//This code produces text files of population mean fitness over time 
//the files for fitness values and initial sequence pool are enclosed in the same folder
// for any inquiry please contact d.morady@gmail.com
//Param 0 iteration (seed)
//Param 1 Mutation rate
//Param 2 recombination rate
//Param 3 population size
//Param 4 Maximum generations (simulation length)
//Param 5 The initial sequence from which the simulation starts
//Param 6 relative weight given to main fitness effects 



public class HIVproject {

    /**
     * @param args
     */
 
    public static void main(String[] args) {
        // TODO Auto-generated method stub
        double [] param=new double[7] ;
 
        param[0]=Integer.parseInt(args[0],10);
        param[1]=Double.parseDouble(args[1]);
        param[2]=Double.parseDouble(args[2]);
        param[3]=Integer.parseInt(args[3],10);
        param[4]=Integer.parseInt(args[4],10);
          param[5]=Integer.parseInt(args[5],10);
          param[6]=Integer.parseInt(args[6],10);
 
        int iteration=1;
        double mutation=0;
 
     
        if(((int) param[1])==6){mutation=0.00001;}
           else if(((int) param[1])==5){mutation=0.00005;}
           else if(((int) param[1])==4){mutation=0.001;}
           else if(((int) param[1])==3){mutation=0.01;}
           else if(((int) param[1])==2){mutation=0.1;}
           else if(((int) param[1])==1){mutation=1;}
        
        
        
        
        
        
        double recombination=0;
        if(((int) param[2])==6){recombination=0;}
        else if(((int) param[2])==5){recombination=0.025;}
        else if(((int) param[2])==4){recombination=0.05;}
        else if(((int) param[2])==3){recombination=0.1;}
        else if(((int) param[2])==2){recombination=0.2;}
        else if(((int) param[2])==1){recombination=0.3;}

        int populationsize=(int) param[3];
        int generation=(int) param[4];
        int seq=(int) param[5];
        double alpha=0;
        if(((int) param[6])==3){alpha=1;}

        PrintWriter meanfitnessf=MakeFile("Meanfitnessseq"+param[5]+"mut"+param[1]+"rec"+param[2]+"alp"+param[6]+"itr"+param[0]);
        Date d = new Date();
        RandomEngine engine = new DRand((int) param[0]);

////////////////////////////////////////

//////////////////////////////////////
        String inn="PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPISPIETVPVKLKPGMDGPKVKQWPSTEEKIKALVEICTEMEKEGKISKIGPENPYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKQKKSVTVLDVGDAYFSVPLDKDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIFQCSMTKILEPFRKQNPDIVIYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGFTTPDKKHQKEPPFLWMGYELHPDKWTVQPIVLPEKDSWTVNDIQKLVGKLNWASQIYAGIKVRQLCKLLRGTKALTEVVPLTEEAELELAE";
 
        ArrayList <Integer> recompos=new ArrayList <Integer> ();
        ArrayList <ArrayList<String>> alleles=new ArrayList <ArrayList<String>> ();

        String input="loci.txt";
        BufferedReader in=null;
        try{
            in=new BufferedReader(new FileReader(input));
        }
        catch(Exception e){System.out.println(e);}
 
 
 
        recompos.add(0);
        //Recombination Matrix
 
 
        while(true){
            try{
                String holder=in.readLine();
                if(holder.matches("null")){break;};
                recompos.add(Integer.parseInt(holder));
            }
            catch(Exception e){
                break;
            }
        }
        ArrayList <String> header=new ArrayList <String> ();
        ArrayList <Integer> headerno=new ArrayList <Integer> ();
 
        String input1="header.txt";
        BufferedReader in1=null;
        try{
            in1=new BufferedReader(new FileReader(input1));
        }
        catch(Exception e){System.out.println(e);}
 

        int counter=0;
        int holder31=0;
        //Recombination Matrix
        String holder1="";
        String holder2="";
        int holder3=0;
          String holder4="";
 
        try{
        holder1=in1.readLine();
        holder2=holder1.substring(0,holder1.length()-2);
        holder4=holder1.substring(holder1.length()-1,holder1.length());
        holder3=Integer.parseInt(holder2.substring(2, holder2.length()));
        }
        catch(Exception e){}

 
 
        String holder11="";
        String holder21="";
        String holder41="";
 
 
        while(true){
            try{
                holder11=in1.readLine();
                if(holder11.matches("null")){break;};
                holder21=holder11.substring(0,holder11.length()-2);
                holder41=holder11.substring(holder11.length()-1,holder11.length());
      
                holder31=Integer.parseInt(holder21.substring(2, holder21.length()));
                if(holder31==holder3){header.add(holder4);counter++;}else{header.add(holder4);alleles.add(header);header=new ArrayList <String> ();headerno.add(counter+1);counter=0;};
                holder3=holder31;
                holder4=holder41;
      
      
            }
            catch(Exception e){
                break;
            }
        }

        ArrayList <Integer> headerid=new ArrayList <Integer> ();
        ArrayList <Integer> headerno1=new ArrayList <Integer> ();
        ArrayList <ArrayList<String>> alleles1=new ArrayList <ArrayList<String>> ();
 
 
 
 
        for(int i=0;i<headerno.size();i++){
            if(headerno.get(i)>1){
                headerid.add(i);
                headerno1.add(headerno.get(i));
                alleles1.add(alleles.get(i));
                for(int j=0;j<alleles.get(i).size();j++){
            //        System.out.print(alleles.get(i).get(j));
                    }
            //    System.out.println(headerno.get(i));
            }
        }

 
 
        double[][] fitnessm=new double [1859][1859];
        String input11="fitness.txt";
        BufferedReader in11=null;
        try{
             in11=new BufferedReader(new FileReader(input11));
            for(int i=0;i<1859;i++){
       
                String[] holder=in11.readLine().split("\t");
                for(int j=0;j<1859;j++){
                    fitnessm[i][j]=Double.parseDouble(holder[j]);
                }
            }
        }
        catch(Exception e){System.out.println(e);}
 
 
        double[][] fitnessint=new double [1859][1859];
        double[][] fitnessall=new double [1][1859];
        String input12="fitnessall.txt";
        BufferedReader in12=null;
        try{
             in12=new BufferedReader(new FileReader(input12));
             for(int i=0;i<1859;i++){
             //    fitnessall[0][i]=0;
                fitnessall[0][i]=Double.parseDouble(in12.readLine());
             }
        }
        catch(Exception e){
         //   System.out.println("didi ridi");
            System.out.println(e);}
 
 
        String input13="fitnessint.txt";
        BufferedReader in13=null;
        try{
             in13=new BufferedReader(new FileReader(input13));
             for(int i=0;i<1859;i++){
                 String[] holder=in13.readLine().split("\t");
                 for(int j=0;j<1859;j++){                  
                         fitnessint[i][j]=Double.parseDouble(holder[j]);
                 }
             }
        }
        catch(Exception e){System.out.println(e);}
        String initialsequence="ISequences.txt";
        String inputt="";
        try{
            in11=new BufferedReader(new FileReader(initialsequence));
        for(int u=0;u<seq;u++){
            inputt=in11.readLine();
        }
 
        }catch(Exception e){System.out.println(e);}
 
    //    System.out.println(inputt);
        String inn1=reversetranslation(headerno,alleles,inputt);
        for(int hh=0;hh<iteration;hh++){
 
 
        String inputbi=reversetranslation(headerno,alleles,inputt);
 
        String[] pop=new String[populationsize];
        double[] fitness=new double[populationsize];
        double nolder=fitnesscalculator(fitnessm,headerno,alleles,inputbi);//fitnesscalculator(fitnessm,headerno,alleles,inputbi);
    //    double nolder=fitnesscalculator2(fitnessall,fitnessint,headerno,alleles, inputbi,alpha);//fitnesscalculator(fitnessm,headerno,alleles,inputbi);
  //     System.out.println(nolder);
 
        for(int i=0; i<populationsize;i++){
            pop[i]=inputbi;
            fitness[i]=nolder;
        }

 
 
 
        for(int h=1;h<generation;h++){
             double Sum=0;
             for(int i=0;i<populationsize;i++){Sum=Sum+fitness[i];};
    
    
            String[] poptemp=new String[populationsize];
            double[] fitnesstemp=new double[populationsize];
        //    sequence.println(h);
      //      fitnessf.println(h);
 
    
   
    
            //int nosample=5;
            double[] matrixy=new double[populationsize];
            double[] matrixx=new double[populationsize];
    
            double[] matrix=new double[populationsize];
            for(int i=0;i<populationsize;i++){
                double tempdouble=Sum*engine.nextDouble();
                matrix[i]=tempdouble;
            }
            Arrays.sort(matrix);
            //System.out.println("sasa");
 
    
            int counter11=0;
            int counter22=0;
    
            double polder=0;
            polder=polder+fitness[counter11];
    
    
    
            while(counter22<matrix.length){
                int sequence1=0;
                int sequence2=0;
                if(polder>=matrix[counter22]){
            
                    sequence1=counter22;
                    matrixy[counter22]=counter11;
                    matrixx[counter22]=counter22;
                    poptemp[counter22]=pop[counter11];
                    fitnesstemp[counter22]=fitness[counter11];
            //        System.out.println("********");
                    counter22++;
                }
                else{
                    counter11++;
                    polder=polder+fitness[counter11];
                }
                        
                if(counter22==(matrix.length)){break;};
        
            }
            double kolder=(populationsize*mutation);
            if(kolder<1){
                double tempo=engine.nextDouble();
                if(tempo<kolder){kolder=1;}
                else{kolder=0;}
            }
            
            
            
            
            int[] matrixm=new int[(int) kolder];
            for(int i=0;i<kolder;i++){
                double tempdouble=populationsize*engine.nextDouble();
                matrixm[i]=(int) tempdouble;
            }
            Arrays.sort(matrixm);
    
            for(int i=0;i<kolder;i++){
                poptemp[matrixm[i]]=mutation(poptemp[matrixm[i]],headerid,recompos,engine);
                fitnesstemp[matrixm[i]]= fitnesscalculator(fitnessm,headerno,alleles,poptemp[matrixm[i]]);
            }
      //      System.out.println(Arrays.toString(matrixm));
    
            kolder=(int) (populationsize*recombination*2);
            int[] matrixr=new int[(int) kolder];
            int hol=0;
            while(hol<kolder){
                int chef=0;
                double tempdouble=populationsize*engine.nextDouble();
                int chofter=(int) tempdouble;
                for(int i=0;i<hol;i++){
                    if(chofter==matrixr[i]){
                        chef=1;
                        break;
                    }
                
                }
                if(chef==0){
                    matrixr[hol]=chofter;
                    hol++;
               
                }
            }
           
    
          //  for(int i=0;i<kolder;i++){
          //      double tempdouble=populationsize*engine.nextDouble();
         //       matrixr[i]=(int) tempdouble;
        //    }
    //        Arrays.sort(matrixr);
           
          
    
            for(int i=0;i<kolder;i=i+2){
                if(!(poptemp[matrixr[i]].equals(poptemp[matrixr[i+1]]))){
                String[] tempd=recombination(poptemp[matrixr[i]], poptemp[matrixr[i+1]],recompos,engine);
                poptemp[i]=tempd[0];
               fitnesstemp[i]=fitnesscalculator(fitnessm,headerno,alleles,tempd[0]);
            //   fitnesstemp[i]=fitnesscalculator2(fitnessall,fitnessint,headerno,alleles, tempd[0],alpha);//fitnesscalculator(fitnessm,headerno,alleles,tempd[0]);
                poptemp[i+1]=tempd[1];
              fitnesstemp[i+1]=fitnesscalculator(fitnessm,headerno,alleles,tempd[1]);
            //    fitnesstemp[i+1]=fitnesscalculator2(fitnessall,fitnessint,headerno,alleles, tempd[1],alpha);//fitnesscalculator(fitnessm,headerno,alleles,tempd[1]);
                }
            }
          //  System.out.println(Arrays.toString(matrixr));
    
           
            pop=poptemp;
           fitness=fitnesstemp;
           meanfitnessf.println(h);
           meanfitnessf.println(Mean(fitness));
    //       System.out.println("");
 
        }
        meanfitnessf.println("111111");
        }
    }
 
    public static double[] bubbleSort(double[] array) {
        boolean swapped = true;
        int j = 0;
        double tmp;
        while (swapped) {
            swapped = false;
            j++;
            for (int i = 0; i < array.length - j; i++) {
                if (array[i] > array[i + 1]) {
                    tmp = array[i];
                    array[i] = array[i + 1];
                    array[i + 1] = tmp;
                    swapped = true;
                }
            }
        }
        return array;
 
    }
 
 
    public static double Mean(double[] input){
        double sum=0;
        for(int i=0;i<input.length;i++){
            sum=sum+input[i];
        }
        return sum/input.length;
    }
 
    public static double fitnesscalculator1(double[][] fitnessm,ArrayList <Integer> headerno, ArrayList <ArrayList<String>> alleles,String inputstring) {
        double fitness=0;
     //  String test=reversetranslation(headerno,alleles,inputstring);
    //    System.out.println(test);
    //    System.out.println(test.length());
        double[][] horizontal=new double[1][1859];
        double[][] vertical=new double[1859][1];
 
        for(int i=0;i<1859;i++){
            double holder=Double.parseDouble(inputstring.substring(i, i+1));
            horizontal[0][i]=holder;
            vertical[i][0]=holder;
      //      System.out.print(holder+" ");
        }
        return fitness;
    }
 
 
    public static double fitnesscalculator2(double[][] fitnessall,double[][] fitnessint,ArrayList <Integer> headerno, ArrayList <ArrayList<String>> alleles,String inputstring, double alpha)  {
        double[][] horizontal=new double[1][1859];
        double[][] vertical=new double[1859][1];
        for(int i=0;i<1859;i++){
            double holder=Double.parseDouble(inputstring.substring(i, i+1));
            horizontal[0][i]=holder;
            vertical[i][0]=holder;
        }
        double res=multipl(fitnessall,vertical)[0][0]+alpha*multiply(multiply(horizontal,fitnessint),vertical)[0][0];
        res=Math.exp(res);
        return res;
    }
 
    public static double variance(double[] v) {
        double mu = Mean(v);
        double sumsq = 0.0;
        for (int i = 0; i < v.length; i++)
          sumsq += sqr(mu - v[i]);
        return sumsq / (v.length);
        // return 1.12; this was done to test a discrepancy with Business
        // Statistics
      }

    public static double sqr(double x) {
        return x*x;
      }
    public static double fitnesscalculator(double[][] fitnessm,ArrayList <Integer> headerno, ArrayList <ArrayList<String>> alleles,String inputstring) {
        double fitness=0;
     //  String test=reversetranslation(headerno,alleles,inputstring);
    //    System.out.println(test);
    //    System.out.println(test.length());
        double[][] horizontal=new double[1][1859];
        double[][] vertical=new double[1859][1];
 
        for(int i=0;i<1859;i++){
            double holder=Double.parseDouble(inputstring.substring(i, i+1));
            horizontal[0][i]=holder;
            vertical[i][0]=holder;
      //      System.out.print(holder+" ");
        }
 
      //  System.out.println(Arrays.toString(multiply(horizontal,fitnessm)[0]));
    //    fitness=multiply(multiply(horizontal,fitnessm),vertical)[0][0]+1;
 
       double[][] vert=new double[1859][1];
    //   vert=multiply(horizontal,fitnessm);
 
 
    //   double sum=0;
     //  for(int i=0;i<1859;i++){
      //      sum=sum+vert[0][i];
       // }
 
 
     //  System.out.println("salam");
      // System.out.println(multiply(horizontal,fitnessm)[0].length);
      // System.out.println(sum);
    //   System.out.println(sum);
       fitness=multiply(multiply(horizontal,fitnessm),vertical)[0][0];
       //System.out.println(multiply(multiply(horizontal,fitnessm),vertical).length);
      // System.out.println(multiply(multiply(horizontal,fitnessm),vertical)[0].length);
       fitness=Math.exp(fitness);
       //System.out.println( fitness);
       //System.out.println("salam");
     //   System.out.println(multiply(multiply(horizontal,fitnessm),vertical)[0][0]);
        return fitness;
    }
    public static void printmatrixd(double[] arg) {
        for(int i=0;i<arg.length;i++){
      //      System.out.println(arg[i]);
        }
    }
    public static void printmatrix(String[] arg) {
        for(int i=0;i<arg.length;i++){
      //      System.out.println(arg[i]);
        }
    }
 
 
    public static double[][] multiply(double a[][], double b[][]){//a[m][n], b[n][p]
           if(a.length == 0) return new double[0][0];
           if(a[0].length != b.length) return null; //invalid dims
 
           int n = a[0].length;
           int m = a.length;
           int p = b[0].length;
 
           double ans[][] = new double[m][p];
 
           for(int i = 0;i < m;i++){
              for(int j = 0;j < p;j++){
                 for(int k = 0;k < n;k++){
                    ans[i][j] += a[i][k] * b[k][j];
                 }
              }
           }
           return ans;
        }
 
 
    public static double[][] multipl(double a[][], double b[][]) {
 
          int aRows = a.length,
              aColumns = a[0].length,
              bRows = b.length,
              bColumns = b[0].length;
          if ( aColumns != bRows ) {throw new IllegalArgumentException("A:Rows: " + aColumns + " did not match B:Columns " + bRows + ".");}
          double[][] resultant = new double[aRows][bColumns];
          for(int i = 0; i < aRows; i++) { // aRow
            for(int j = 0; j < bColumns; j++) { // bColumn
              for(int k = 0; k < aColumns; k++) { // aColumn
                resultant[i][j] += a[i][k]*b[k][j];
    
              }
     
            }
          }
     //     System.out.println(resultant[0][0]);
          return resultant;
        }
    public static String[] recombination(String s1, String s2,ArrayList <Integer> recompos, RandomEngine engine){
        String[] array=new String[2];
        String s1t="";
        String s2t="";
           int holder=(int) (engine.nextDouble()*recompos.size()-1);
           s1t=s1.substring(0, recompos.get(holder+1)).concat(s2.substring(recompos.get(holder+1), s1.length()));
           s2t=s2.substring(0, recompos.get(holder+1)).concat(s1.substring(recompos.get(holder+1), s2.length()));
           array[0]=s1t;
           array[1]=s2t;
           return array;
          }
    public static String mutation(String inputtrans, ArrayList <Integer> headerid,ArrayList <Integer> recompos, RandomEngine engine){
 
 
            int mut= (int) (engine.nextDouble()*headerid.size());
           String tempmut=inputtrans.substring(recompos.get(headerid.get(mut)), recompos.get(headerid.get(mut)+1));
           ArrayList <Integer> arrtemp=new ArrayList <Integer> ();
           int polder=0;
 
           for(int j=0;j<tempmut.length();j++){
               String tempm=tempmut.substring(j, j+1);
               if(!tempm.equals("1")){arrtemp.add(j);}else{polder=j;};
           }
 
           int temp1=(int) (engine.nextDouble()*arrtemp.size());
         //  System.out.println(temp1);
           StringBuffer temp2=new StringBuffer(tempmut);
           temp2=temp2.replace(arrtemp.get(temp1), arrtemp.get(temp1)+1, "1");
           temp2=temp2.replace(polder, polder+1, "0");
        //   System.out.println(tempmut);
        //   System.out.println(temp2);
 
           StringBuffer outtemp=new StringBuffer(inputtrans);
           outtemp=outtemp.replace(recompos.get(headerid.get(mut)), recompos.get(headerid.get(mut)+1), temp2.toString());
        //   System.out.println(outtemp);
 
 
           return outtemp.toString();
 
    }
     public static String reversetranslation(ArrayList <Integer> headerno, ArrayList <ArrayList<String>> alleles, String inputstring) {
        String inputtrans="";
 
        for(int i=0;i<headerno.size();i++){
            int[] temp=new int[headerno.get(i)];
            String temp2=inputstring.substring(i, i+1);
            String temp1="";
            for(int j=0;j<headerno.get(i);j++){
                if(alleles.get(i).get(j).equals(temp2)){temp1=temp1.concat("1");}
                else{temp1=temp1.concat(Integer.toString(temp[j]));}
            }
            inputtrans=inputtrans.concat(temp1);
        }
        return inputtrans;
 
    }
     public static PrintWriter MakeFile(String input) {
            File file = new File(input);
            file.delete();
            File file1 = new File(input);
            FileWriter fw=null;
            try {
                fw = new FileWriter(file1, true);
            } catch (IOException e) {
                // TODO Auto-generated catch block
           //     System.out.println(e);
                }
            BufferedWriter bw = new BufferedWriter(fw);
            PrintWriter out = new PrintWriter(bw, true);
            return out;
        }
     public static String translation(ArrayList <Integer> recompos, ArrayList <ArrayList<String>> alleles, String inputtrans) {
         String reconstruct="";
         for(int i=0;i<recompos.size()-1;i++){
             String temp=inputtrans.substring(recompos.get(i), recompos.get(i+1));
             for(int j=0;j<temp.toCharArray().length;j++){
                 char hold=temp.toCharArray()[j];
                 if(Character.toString(hold).equals("1")){
                     reconstruct=reconstruct.concat(alleles.get(i).get(j));
                 }
             }
         }
         return reconstruct;
     }
}