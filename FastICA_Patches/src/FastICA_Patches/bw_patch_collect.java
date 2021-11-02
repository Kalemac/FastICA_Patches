package FastICA_Patches;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.Scanner;
import java.util.concurrent.ThreadLocalRandom;
import javax.imageio.ImageIO;

import org.ejml.data.Matrix;
import org.fastica.BelowEVFilter;
import org.fastica.CompositeEVFilter;
import org.fastica.FastICA;
import org.fastica.FastICAConfig;
import org.fastica.ProgressListener;
import org.fastica.FastICAConfig.Approach;
import org.fastica.FirstEVFilter;
import org.fastica.GaussCFunction;
import org.fastica.PercentageEVFilter;
import org.fastica.Power3CFunction;
import org.fastica.ProgressListener.ComputationState;
import org.fastica.SortingEVFilter;
import org.fastica.TanhCFunction;

import com.opencsv.CSVReader;

//import org.fastica.math.Matrix;

import com.opencsv.CSVWriter;
import java.io.FileWriter;

import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;


public class bw_patch_collect {

	static double findMatrixMean(double m[][][], int w, int h, int d) {

		double sum = 0;
		for(int i = 0; i < w; i++) {
			for(int j = 0; j < h; j++) {
				for(int k = 0; k < d; k++) {
					sum += m[i][j][k];
				}
			}
		}
		return sum / (w * h * d);
	}
	
	static double findMatrixMean(double m[][], int w, int h) {

		double sum = 0;
		for(int i = 0; i < w; i++) {
			for(int j = 0; j < h; j++) {
				sum += m[i][j];
			}
		}
		return sum / (w * h);
	}

	static double findArrayMean(double array[]) {
		double sum = 0;

		for(int i = 0; i < array.length; i++)
			sum += array[i];

		return sum / array.length;
	}

	static double findVariance(double matrix[][][], int w, int h, int d, double mean) {
		double sum = 0;
		
		double[][][] m = new double[matrix.length][matrix[0].length][];
		for(int i = 0; i < matrix.length; i++) {
			for(int j = 0; j < matrix[0].length; j++) {
				m[i][j] = matrix[i][j].clone();
			}
		}
		
		for(int i = 0; i < w; i++) {
			for(int j = 0; j < h; j++) {
				for(int k = 0; k < d; k++) {
					m[i][j][k] -= mean;
					m[i][j][k] *= m[i][j][k];
					sum += m[i][j][k];
				}
			}
		}
			
		return sum / (w * h * d);
		
	}
	
	static double findVariance(double matrix[][], int w, int h, double mean) {
		double sum = 0;

		double [][] m = new double[matrix.length][];
		for(int i = 0; i < matrix.length; i++)
			m[i] = matrix[i].clone();

		for(int i = 0; i < w; i++) {
			for(int j = 0; j < h; j++) {
				m[i][j] -= mean;
				m[i][j] *= m[i][j];
			}
		}

		for(int i = 0; i < w; i++) {
			for(int j = 0; j < h; j++) {
				sum += m[i][j];
			}
		}

		return sum / (w * h);
	}

	static double[][] subMatrix(double[][] matrix, int row, int col, int width){
		double[][] sub = new double[width][width];

		for(int i = row; i < row+width; i++) {
			for(int j = col; j < col+width; j++) {
				int x = i - row;
				int y = j - col;

				sub[x][y] = matrix[i][j];
			}
		}

		return sub;
	}
	
	static double[][][] subMatrix(double[][][] matrix, int row, int col, int width){
		double [][][] sub = new double[width][width][3];
		
		for(int i = row; i < row+width; i++) {
			for(int j = col; j < col+width; j++) {
				int x = i - row;
				int y = j - col;
				for(int k = 0; k < 3; k++) {
					sub[x][y][k] = matrix[i][j][k];
				}
			}
		}
		
		return sub;
	}
	
	static double[] reshapeMatrix(double[][] matrix, int size) {
		double[] reshaped = new double[size];
		int index = 0;

		for(int i = 0; i < matrix.length; i++) {
			for(int j = 0; j < matrix[i].length; j++) {
				reshaped[index] = matrix[i][j];
				index++;
			}
		}

		return reshaped;
	}
	
	private static double[] reshapeMatrix(double[][][] matrix, int size) {
		double[] reshaped = new double[size];
		int index = 0;

		for(int i = 0; i < matrix.length; i++) {
			for(int j = 0; j < matrix[i].length; j++) {
				for(int k = 0; k < 3; k ++) {
					reshaped[index] = matrix[i][j][k];
					index++;
				}
			}
		}

		return reshaped;
	}
	
	static double[][] reshapeArray(double array[], int row, int col){
		double[][] rMatrix = new double[row][col];
		
		int i = 0;
		
		for (int r = 0; r < row; r++) {
			for (int c = 0; c < col; c++) {
				rMatrix[r][c] = array[i++];
			}
		}
		
		return rMatrix;
	}

	static double findMax(double[] array) {
		double max = array[0];
		for(int i = 1; i < array.length; i++) {
			if(array[i] > max)
				max = array[i];
		}

		return max;
	}

	static double findMax(double[][] array) {
		double max = array[0][0];
		for(int i = 0; i < array.length; i++) {
			for(int j = 0; j < array[i].length; j++) {
				if(array[i][j] > max)
					max = array[i][j];
			}
		}

		return max;
	}

	static double findMin(double[] array) {
		double min = array[0];
		for(int i = 1; i < array.length; i++) {
			if(array[i] < min)
				min = array[i];
		}

		return min;
	}
	
	static double findMin(double[][] array) {
		double min = array[0][0];
		for(int i = 0; i < array.length; i++) {
			for(int j = 0; j < array[i].length; j++) {
				if(array[i][j] < min)
					min = array[i][j];
			}
		}

		return min;
	}
	
	static void print2dArray(double[][] array) {
		for(int i = 0; i < array.length; i++) {
			for(int j = 0; j < array[i].length; j++) {				
				System.out.print(new DecimalFormat("#.###").format(array[i][j]) + " ");
			}
			System.out.println();
		}
	}
	
	static void print3dArray(double[][][] array) {
		for(int i = 0; i < array.length; i++) {
			for(int j = 0; j < array[i].length; j++) {
				for(int k = 0; k < 3; k++) {
					System.out.print(new DecimalFormat("#.###").format(array[i][j][k]) + " ");
				}
				System.out.print("//");
			}
			System.out.println();
		}
	}
	
	public static double[][] getPatches(int nPatches, int pWidth, String filePath){

		int numPatches = nPatches;
        int patchWidth = pWidth;


        int numPixels = patchWidth * patchWidth;
        int patchCount = 0;
        int tryCount = 0;

        double[][] patchSample = new double[patchWidth][patchWidth];
        double[] patch = new double[numPixels];
        double[][] imgPatches = new double[numPixels][numPatches];

        //Bitmap image = bitmap;

        //Load image to be processed
        BufferedImage image = null;
        try {
            image = ImageIO.read(bw_patch_collect.class.getResource(filePath));
        } catch (IOException e) {
            e.printStackTrace();
        }

        int imageHeight = image.getHeight();
        int imageWidth = image.getWidth();

        double[][] imageArray = new double[imageWidth][imageHeight];	//create pixel array

        //Convert individual color pixel to greyscale
        for(int i = 0; i < imageWidth; i++) {
            for(int j = 0; j < imageHeight; j++) {
                int rgb = image.getRGB(i,j);
                int r = (rgb >> 16) & 0xFF;
                int g = (rgb >> 8) & 0xFF;
                int b = (rgb & 0xFF);

                //THIS IS A LINEAR APROXIMATION
                imageArray[i][j] = Math.round(0.3*r + .59*g + .11*b);
                //UPDATE THIS LATER FOR TRUE GREYSCALE
            }
        }

        //calculate mean of matrix
        double mean = findMatrixMean(imageArray, imageWidth, imageHeight);

        //subtract mean from each value in the matrix
        for(int i = 0; i < imageWidth; i++) {
            for(int j = 0; j < imageHeight; j++) {
                imageArray[i][j] -= mean;
            }
        }
        //calculate Variance of new matrix
        double var = findVariance(imageArray, imageWidth, imageHeight, findMatrixMean(imageArray, imageWidth, imageHeight));
        //calculate standard deviation of new matrix
        double std = Math.sqrt(var);

        //divide each value in the matrix by the standard deviation
        for(int i = 0; i < imageWidth; i++) {
            for(int j = 0; j < imageHeight; j++) {
                imageArray[i][j] /= std;
            }
        }

        while(patchCount < numPatches && tryCount < numPatches) {
            tryCount += 1;

            //Create random X and Y value for patch
            int px = new Random().nextInt(imageWidth - patchWidth + 1);
            int py = new Random().nextInt(imageHeight - patchWidth + 1);

//          int px = ThreadLocalRandom.current().nextInt(0, imageWidth - patchWidth + 1);
//          int py = ThreadLocalRandom.current().nextInt(0, imageHeight - patchWidth + 1);

            //px = 200;
            //py = 150;

            patchSample = subMatrix(imageArray, px, py, patchWidth);	//grab patch from image array
            double patchStd = Math.sqrt(findVariance(patchSample, patchWidth, patchWidth, findMatrixMean(patchSample, patchWidth, patchWidth)));	//find standard deviation of the patch

            if(patchStd > 0.0) {
                patch = reshapeMatrix(patchSample, numPixels);	//Reshape 2-dimensional array to 1-dimensional array
                double patchMean = findArrayMean(patch);		//calculate mean of the array
                for(int i = 0; i < patch.length; i++) {
                    patch[i] = patch[i] - patchMean;			//subtract mean from each element in array
                }
                for(int i = 0; i < patch.length; i++) {
                    imgPatches[i][patchCount] = patch[i];		//copy individual patch into array of patches
                }
                patchCount += 1;
            }
        }

        return imgPatches;
    }
	
	public static double[][] getPatchesColor(int nPatches, int pWidth, String filePath){
		int numPatches = nPatches;
		int patchWidth = pWidth;


		int numPixels = patchWidth * patchWidth;
		int patchCount = 0;
		int tryCount = 0;

		double[][][] patchSample = new double[patchWidth][patchWidth][3];
		//double[] patch = new double[numPixels];
		double[][] imgPatches = new double[numPixels * 3][numPatches];
		
		//Load image to be processed
		BufferedImage image = null;
		try {
			image = ImageIO.read(bw_patch_collect.class.getResource(filePath));
		} catch (IOException e) {
			e.printStackTrace();
		}

		int imageHeight = image.getHeight();
		int imageWidth = image.getWidth();
		
		double[][][] imageArray = new double[imageWidth][imageHeight][3];
		
		for(int i = 0; i < imageWidth; i++) {
			for(int j = 0; j < imageHeight; j++) {
				int rgb = image.getRGB(i,j);
				imageArray[i][j][0] = (rgb >> 16) & 0xFF;	//red
				imageArray[i][j][1] = (rgb >> 8) & 0xFF;	//green
				imageArray[i][j][2] = (rgb & 0xFF);			//blue

			}
		}
		
		//calculate mean of matrix
		double mean = findMatrixMean(imageArray, imageWidth, imageHeight, 3);
		
		//subtract mean from each value in the matrix
		for(int i = 0; i < imageWidth; i++) {
			for(int j = 0; j < imageHeight; j++) {
				for(int k = 0; k < 3; k++) {
					imageArray[i][j][k] -= mean;
				} 
			}
		}
		
		//calculate Variance of new matrix
		double var = findVariance(imageArray, imageWidth, imageHeight, 3, findMatrixMean(imageArray, imageWidth, imageHeight, 3));
		//calculate standard deviation of new matrix
		double std = Math.sqrt(var);
		
		for(int i = 0; i < imageWidth; i++) {
			for(int j = 0; j < imageHeight; j++) {
				for(int k = 0; k < 3; k++) {
					imageArray[i][j][k] /= std; 
				}
			}
		}
		
		while(patchCount < numPatches && tryCount < numPatches) {
			tryCount += 1;
			
			int px = ThreadLocalRandom.current().nextInt(0, imageWidth - patchWidth + 1);		
			int py = ThreadLocalRandom.current().nextInt(0, imageHeight - patchWidth + 1);
			
//			int px = 1;
//			int py = 1;
			
			patchSample = subMatrix(imageArray, px, py, patchWidth);
			double patchStd = Math.sqrt(findVariance(patchSample, patchWidth, patchWidth, 3, findMatrixMean(patchSample, patchWidth, patchWidth, 3)));
			
			if(patchStd > 0.0) {
				double[] patch = reshapeMatrix(patchSample, numPixels * 3);	//Reshape 2-dimensional array to 1-dimensional array
				//double[] patch = Arrays.stream(patchSample).flatMap(Arrays::stream).flatMapToDouble(Arrays::stream).toArray();
				//System.out.println(Arrays.toString(patch));
				double patchMean = findArrayMean(patch);			//calculate mean of the array
				for(int i = 0; i < patch.length; i++) {
					patch[i] = patch[i] - patchMean;				//subtract mean from each element in array
				}
				
				for(int i = 0; i < patch.length; i++) {
					imgPatches[i][patchCount] = patch[i];		//copy individual patch into array of patches
				}
				patchCount += 1;
			}
			
		}
		
		return imgPatches;
	}
	

	public static void showPatches(double[][] prePatches, int showPatchNum, boolean display, String name) {
		
		double[][] patches = prePatches;
		
		int dataDim = patches.length;				//number of patches in the 2-dimensional array
		int patchWidth = (int)Math.sqrt(dataDim);	//calculate  width of a patch

		double[][] displayPatch = new double[dataDim][showPatchNum];	//final 2-dimensional with all patches and gaps between patches

		double[] patch = new double[dataDim];

		for(int i = 0; i < showPatchNum; i++) {
			int patch_i = i;
			for(int j = 0; j < dataDim; j++) {
				patch[j] = patches[j][patch_i];	
			}
			double pmax = findMax(patch);
			double pmin = findMin(patch);
			
			//Not entirely sure why we need this if statement since it should ALWAYS run
			//But I will keep it in because the base code does
			if(pmax > pmin) {	
				for(int k = 0; k < patch.length; k++) {
					patch[k] = (patch[k] - pmin) / (pmax - pmin);	//convert values so they range from 0 to 1
				}
			}

			for(int j = 0; j < dataDim; j++) {
				displayPatch[j][i] = patch[j];
			}
		}

		int bw = 5;			//Border Width
		int pw = patchWidth;

		int patchesY = (int)Math.sqrt(showPatchNum);
		int patchesX = (int)Math.ceil(((double)showPatchNum) / patchesY);

		double[][] patchImg = new double[(pw + bw) * patchesX - bw][patchesY * (pw + bw) - bw];

		double dpMax = findMax(displayPatch);

		for(int i = 0; i < patchImg.length; i++) {
			for(int j = 0; j < patchImg[i].length; j++) {
				patchImg[i][j] = dpMax;
			}
		}

		for (int i = 0; i < showPatchNum; i++) {
			int y_i = Math.floorDiv(i,  patchesY);
			int x_i = i % patchesY;
			
			double[][] reshaped = new double[pw][pw];
			int x = 0;
			for (int j = 0; j < reshaped.length; j++) {
				for (int j2 = 0; j2 < reshaped[j].length; j2++) {
					reshaped[j][j2] = displayPatch[x++][i];
				}				
			}
			
			int rX = 0;
			for(int row = x_i * (pw + bw); row < x_i * (pw + bw) + pw; row++) {
				int rY = 0;
				for(int col = y_i * (pw + bw); col < y_i * (pw + bw) + pw; col++) {
					patchImg[row][col] = reshaped[rX][rY];
					rY++;					
				}
				rX++;
			}
		}
			
		
		if(display) {
			int xLength = patchImg.length;
			int yLength = patchImg[0].length;
			BufferedImage buff = new BufferedImage(xLength, yLength, BufferedImage.TYPE_3BYTE_BGR);		//create blank image to be written to
			
			for(int x = 0; x < xLength; x++) {
			    for(int y = 0; y < yLength; y++) {
			    	int value = (int)(patchImg[x][y] * 255);	//convert pixel value to 0 to 255 range
			        int rgb = value << 16 | value << 8 | value;	//assign red, green, and blue values to the same to create a greyscale image
			        buff.setRGB(x, y, rgb);		//assign pixel to rbg value
			    }
			}
			try {
				ImageIO.write(buff, "png", new File(String.valueOf(name) + ".png"));	//save image to project directory
			} catch (IOException e) {
				e.printStackTrace();
			}
			
		}
	}
	
	public static void showPatchesColor(double[][] prePatches, int showPatchNum, boolean display, String name) {
		double[][] patches = prePatches;
		
		int dataDim = patches.length;				//number of patches in the 2-dimensional array
		int patchWidth = (int)Math.sqrt(dataDim/3);	//calculate  width of a patch
		
		double[][] displayPatch = new double[dataDim][showPatchNum];	//final 2-dimensional with all patches and gaps between patches
		
		double[] patch = new double[dataDim];
		
		for (int i = 0; i < showPatchNum; i++) {
			for(int j = 0; j < dataDim; j++) {
				patch[j] = patches[j][i];
			}
			double pmax = findMax(patch);
			double pmin = findMin(patch);
			
			//Not entirely sure why we need this if statement since it should ALWAYS run
			//But I will keep it in because the base code does
			if(pmax > pmin) {	
				for(int k = 0; k < patch.length; k++) {
					patch[k] = (patch[k] - pmin) / (pmax - pmin);	//convert values so they range from 0 to 1
				}
			}
			
			for(int j = 0; j < dataDim; j++) {
				displayPatch[j][i] = patch[j];
			}
			
		}
		
		int bw = 5;
		int pw = patchWidth;
		
		int patchesY = (int)Math.sqrt(showPatchNum);
		int patchesX = (int)Math.ceil(((double)showPatchNum) / patchesY);
		
		double[][][] patchImg = new double[(pw + bw) * patchesX - bw][patchesY * (pw + bw) - bw][3];
		
		//double dpMax = findMax(displayPatch);
		
		for(int i = 0; i < patchImg.length; i++) {
			for(int j = 0; j < patchImg[i].length; j++) {
				for(int k = 0; k < 3; k++) {
					patchImg[i][j][k] = 1.0; 
				}
			}
		}

		
		for (int i = 0; i < showPatchNum; i++) {
			int y_i = Math.floorDiv(i,  patchesY);
			int x_i = i % patchesY;
			
			double[][][] reshaped = new double[pw][pw][3];
			int x = 0;
			
			for (int j = 0; j < reshaped.length; j++) {
				for (int j2 = 0; j2 < reshaped[j].length; j2++) {
					for (int j3 = 0; j3 < 3; j3++) {
						reshaped[j][j2][j3] = displayPatch[x++][i];
					}
				}				
			}
			
			int rX = 0;
			for(int row = x_i * (pw + bw); row < x_i * (pw + bw) + pw; row++) {
				int rY = 0;
				for(int col = y_i * (pw + bw); col < y_i * (pw + bw) + pw; col++) {
					for(int rgb = 0; rgb < 3; rgb++) {
						patchImg[row][col][rgb] = reshaped[rX][rY][rgb];
					}
					rY++;					
				}
				rX++;
			}	
			
		}
		
		if(display) {
			int xLength = patchImg.length;
			int yLength = patchImg[0].length;
			BufferedImage buff = new BufferedImage(xLength, yLength, BufferedImage.TYPE_3BYTE_BGR);		//create blank image to be written to
			
			for(int x = 0; x < xLength; x++) {
			    for(int y = 0; y < yLength; y++) {
			        int rgb = (int)(patchImg[x][y][0] * 255) << 16 | (int)(patchImg[x][y][1] * 255) << 8 | (int)(patchImg[x][y][2] * 255);	//assign red, green, and blue values to the same to create a greyscale image
			        buff.setRGB(x, y, rgb);		//assign pixel to rbg value
			    }
			}
			try {
				ImageIO.write(buff, "png", new File(String.valueOf(name) + ".png"));	//save image to project directory
			} catch (IOException e) {
				e.printStackTrace();
			}
			
		}
		
	}
	
	public static double[][] transpose(double[][] matrix) {
        int m = matrix.length;
        int n = matrix[0].length;
        double[][] transposed = new double[n][m];
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                transposed[i][j] = matrix[j][i];
            }
        }
        return transposed;
    }
	
	
	
	public static void main(String[] args) throws Exception {


		int numPatches = 10000;
		int patchSize = 8;
		int numComponents = 25;		//must be a perfect square ie. 4, 9, 16, 25, etc...
		
		//Let the user know the program has started as well as start the clock to calculate runtime
		System.out.println("Program has Started");
		long startTime = System.nanoTime();
		
		//Get patches from the supplied image
		double [][] prePatches = getPatches(numPatches, patchSize, "/images/redpanda.jpg");
		
		BufferedWriter outputWriter = null;
		  outputWriter = new BufferedWriter(new FileWriter("Pre_SCV.csv"));
		  for (int i = 0; i < prePatches.length; i++) {
			  for (int j = 0; j < prePatches[i].length; j++) {
				  outputWriter.write(Double.toString(prePatches[i][j])+", ");
			  }
		    outputWriter.newLine();
		  }
		  outputWriter.flush();  
		  outputWriter.close();  
		
		
		
		
		//create file name with descriptors
		String prePatch = "pre-patch" + String.valueOf(numPatches) + "_" + String.valueOf(patchSize) + "_" + String.valueOf(numComponents);
		
		//Display the patches before being passed to fastICA
		showPatches(prePatches, numComponents, true, prePatch);
		
		//Create fastICA object then fit the patches
		
//		FastICAConfig config = new FastICAConfig(numComponents, Approach.DEFLATION, 1.0, 1.0e-15, 5000, null);
		
//		CompositeEVFilter filter = new CompositeEVFilter();
	
		//filter.add(new SortingEVFilter(true, true));
		//filter.add(new FirstEVFilter(numComponents));
		
//		filter.add(new BelowEVFilter(0.1, false));
		//filter.add(new FirstEVFilter(numComponents));
		//filter.add(new SortingEVFilter(true, true));
		
//		System.out.println("FastICA is starting...");
		
//		FastICA test = new FastICA(prePatches, config, new TanhCFunction(1.0), filter, new ProgressListener() {
//			public void progressMade( 
//					ComputationState state, 
//					int component, 
//					int iteration, 
//					int maxComps )
//				{
//				}
//			
//		});
		

		//FastICA fastICA = new FastICA();	
		//fastICA.fit(transpose(prePatches), numComponents);
		
		//get the releulting mixing matrix
//		System.out.println("Getting Mixing Matrix...");
//		double[][] mixing = test.getMixingMatrix();
//		double[][] compVec = test.getICVectors();
//		double[][] unMixing = test.getSeparatingMatrix();
		
		//set name for the post-processed patches
		String postPatch = "post-patch" + String.valueOf(numPatches) + "_" + String.valueOf(patchSize) + "_" + String.valueOf(numComponents);
		
//		System.out.println("Mixing: " + mixing.length + ", " + mixing[0].length);
//		System.out.println("Component: " + compVec.length + ", " + compVec[0].length);
//		System.out.println("Un-Mixing: " + unMixing.length + ", " + unMixing[0].length);
//		
//		showPatches(mixing, numComponents, true, "Mixing_Matrix_Test");
//		showPatches(compVec, numComponents, true, "Component_Test");
//		showPatches(transpose(unMixing), numComponents, true, "Un-Mixing_Test");
		
//		Scanner sc = new Scanner(new BufferedReader(new FileReader("filters.csv")));
//		  int rows = 64;
//		  int columns = 36;
//		  double [][] myArray = new double[rows][columns];
//		  while(sc.hasNextLine()) {
//		     for (int i=0; i<myArray.length; i++) {
//		        String[] line = sc.nextLine().trim().split(",");
//		        for (int j=0; j<line.length; j++) {
//		           myArray[i][j] = Double.parseDouble(line[j]);
//		        }
//		     }
//		  }
		
		//display the patches after being passed into fastICA
//		showPatches(myArray, 16, true, "test");
		
		//End clock and print resulting runtime
		long endTime   = System.nanoTime();
		long totalTime = endTime - startTime;
		float sec =  totalTime / 1000000000F;
		System.out.println("Runtime: " + sec + " seconds");
		
		boolean print = false;

		if(print) {
			print2dArray(prePatches);
		}

	}
}
