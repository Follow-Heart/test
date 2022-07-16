//#include "..\\Header\\sim.h"
#include "sim.h"
#include <math.h>

using namespace std;

int pcounter_n0;
int pcounter_c0;
int pcounter_k0;
int pcounter_m0;

int pcounter_tn0;
int pcounter_tc0;
int pcounter_tk0;
int pcounter_tm0;

int pcounter_m1;
int pcounter_c1;
int pcounter_n1;
int pcounter_k1;

int pcounter_tm1;
int pcounter_tc1;
int pcounter_tn1;
int pcounter_tk1;

void initSim()
{
	pcounter_n0 = 0;
	pcounter_c0 = 0;
	pcounter_k0 = 0;
	pcounter_m0 = 0;

	pcounter_tn0 = 0;
	pcounter_tc0 = 0;
	pcounter_tk0 = 0;
	pcounter_tm0 = 0;

	pcounter_m1 = 0;
	pcounter_c1 = 0;
	pcounter_n1 = 0;
	pcounter_k1 = 0;

	pcounter_tm1 = 0;
	pcounter_tc1 = 0;
	pcounter_tn1 = 0;
	pcounter_tk1 = 0;
}

int countDramAccess(std::vector<std::vector<double>> &MAT)
{
	int count = 0;
	for (int i = 0; i < MAT.size(); i++)
	{
		for (int j = 0; j < MAT[i].size(); j++)
		{
			if (MAT[i][j] != 0)
				count++;
		}
	}
	return count;
}

void sim()
{

	char sf[] = "**********************************************************************************************************************************\n";
	char end0[] = "\n***********************************************************************************************************************************\n   ----------------------------------------------------------Sim has done-------------------------------------------------------\n";
	char end1[] = "   ------------------------------------------------------  cycle number is : ";
	char end2[] = "  ---------------------------------------------   \n";
	char end3[] = "   ---------------------------------------------  ";
	char end4[] = "  -------------------------------------------   \n";
	char end5[] = "***********************************************************************************************************************************\n";

	long long cycle = 0;
	long long num_dram_access = 0;
	long long num_sram_access = 0;
	long long prev_cycle = 0;
	long long prev_num_dram_access = 0;
	long long prev_num_sram_access = 0;
	long long compute_latency = 0;
	long long dram_access_latency = 0;

	initDram();
	//assignDram(DATASET);
	initSram();
	initReg();
	initSim();
	configureWorkload();

	int prev_pcounter_n0 = -1;
	int prev_pcounter_c0 = -1;
	int prev_pcounter_k0 = -1;
	int prev_pcounter_m0 = -1;

	int prev_pcounter_tn0 = -1;
	int prev_pcounter_tc0 = -1;
	int prev_pcounter_tk0 = -1;
	int prev_pcounter_tm0 = -1;

	int prev_pcounter_m1 = -1;
	int prev_pcounter_c1 = -1;
	int prev_pcounter_n1 = -1;
	int prev_pcounter_k1 = -1;

	int prev_pcounter_tm1 = -1;
	int prev_pcounter_tc1 = -1;
	int prev_pcounter_tn1 = -1;
	int prev_pcounter_tk1 = -1;

	clock_t begin_sim, end_sim;
	begin_sim = clock();

	if (ACCELERATOR == "GCNAX")
	{
		if (!LOOP_FUSION)
		{
			prev_pcounter_n0 = -1;
			prev_pcounter_c0 = -1;
			prev_pcounter_k0 = -1;
			prev_pcounter_tn0 = -1;
			prev_pcounter_tc0 = -1;
			prev_pcounter_tk0 = -1;
			prev_pcounter_m1 = -1;
			prev_pcounter_c1 = -1;
			prev_pcounter_n1 = -1;
			prev_pcounter_tm1 = -1;
			prev_pcounter_tc1 = -1;
			prev_pcounter_tn1 = -1;

			//off-chip communication: SPMM1, B=X*W
			//load matrix X, W, B to SUBMAT_X, SUBMAT_W, SUBMAT_BO

			for (pcounter_n0 = 0; pcounter_n0 < DIM_N; pcounter_n0 += TILE_SIZE_TN0)
			{
				for (pcounter_c0 = 0; pcounter_c0 < DIM_C; pcounter_c0 += TILE_SIZE_TC0)
				{
					for (pcounter_k0 = 0; pcounter_k0 < DIM_K; pcounter_k0 += TILE_SIZE_TK0)
					{
						prev_num_dram_access = num_dram_access; //record the current number of dram accesses

						if (pcounter_n0 != prev_pcounter_n0 || pcounter_k0 != prev_pcounter_k0)
						{
							//  load SUBMAT_X
							for (int i = 0; i < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); i++)
							{
								SUBMAT_X[i].assign(MAT_X[i + pcounter_n0].begin() + pcounter_k0, MAT_X[i + pcounter_n0].begin() + std::min(pcounter_k0 + TILE_SIZE_TK0, DIM_K));
								//SUBMAT_X[i].resize(TILE_SIZE_TK);
							}
							//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TK, DIM_K - pcounter_k);
							num_dram_access += countDramAccess(SUBMAT_X);
						}

						if (pcounter_k0 != prev_pcounter_k0 || pcounter_c0 != prev_pcounter_c0)
						{
							//  load SUBMAT_W
							for (int i = 0; i < std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0); i++)
							{
								SUBMAT_W[i].assign(MAT_W[i + pcounter_k0].begin() + pcounter_c0, MAT_W[i + pcounter_k0].begin() + std::min(pcounter_c0 + TILE_SIZE_TC0, DIM_C));
								//SUBMAT_W[i].resize(TILE_SIZE_TC0);
							}
							//num_dram_access += std::min(TILE_SIZE_TK, DIM_K - pcounter_k) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							num_dram_access += countDramAccess(SUBMAT_W);
						}

						if (pcounter_n0 != prev_pcounter_n0 || pcounter_c0 != prev_pcounter_c0)
						{
							//  load SUBMAT_BO
							// if SUBMAT_BO is valid, which means it has to write back to DRAM
							if (valid_submat_bo && prev_pcounter_n0 != -1 && prev_pcounter_c0 != -1)
							{
								for (int i = 0; i < std::min(TILE_SIZE_TN0, DIM_N - prev_pcounter_n0); i++)
								{
									for (int j = 0; j < std::min(TILE_SIZE_TC0, DIM_C - prev_pcounter_c0); j++)
									{
										MAT_B[i + prev_pcounter_n0][j + prev_pcounter_c0] = SUBMAT_BO[i][j];
									}
								}
								//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
								num_dram_access += countDramAccess(SUBMAT_BO);
								valid_submat_bo = false;
							}

							//  load SUBMAT_BO
							for (int i = 0; i < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); i++)
							{
								SUBMAT_BO[i].assign(MAT_B[i + pcounter_n0].begin() + pcounter_c0, MAT_B[i + pcounter_n0].begin() + std::min(pcounter_c0 + TILE_SIZE_TC0, DIM_C));
								//SUBMAT_BO[i].resize(TILE_SIZE_TC0);
							}
							//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							num_dram_access += countDramAccess(SUBMAT_BO);
						}

						dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access - prev_num_dram_access) / (DRAM_BANDWIDTH / CLOCK_RATE); //calculate the dram access latency

						//on-chip computation
						// loop order tn0->tk->tc0
						prev_cycle = cycle;
						for (pcounter_tn0 = 0; pcounter_tn0 < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); pcounter_tn0++)
						{
							for (pcounter_tk0 = 0; pcounter_tk0 < std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0); pcounter_tk0++)
							{
								if (SUBMAT_X[pcounter_tn0][pcounter_tk0] != 0)
								{
									if (pcounter_tn0 != prev_pcounter_tn0 || pcounter_tk0 != prev_pcounter_tk0)
									{
										REG_SCALAR = SUBMAT_X[pcounter_tn0][pcounter_tk0];
										num_sram_access += 1;
									}

									for (pcounter_tc0 = 0; pcounter_tc0 < std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0); pcounter_tc0 += MULTIPLIER_NUM)
									{
										// load element of W to registers
										if (pcounter_tk0 != prev_pcounter_tk0 || pcounter_tc0 != prev_pcounter_tc0)
										{
											for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc0); i++)
											{
												REG_VEC1[i] = SUBMAT_W[pcounter_tk0][pcounter_tc0 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}

										// load elements of B to registers
										if (pcounter_tn0 != prev_pcounter_tn0 || pcounter_tc0 != prev_pcounter_tc0)
										{
											// if replaced, the results in registers need write back
											if (valid_vec2 && prev_pcounter_tn0 != -1 && prev_pcounter_tc0 != -1)
											{
												for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - prev_pcounter_tc0); i++)
												{
													SUBMAT_BO[prev_pcounter_tn0][prev_pcounter_tc0 + i] = REG_VEC2[i];
												}
												valid_vec2 = false;
												num_sram_access += MULTIPLIER_NUM;
											}

											for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc0); i++)
											{
												REG_VEC2[i] = SUBMAT_BO[pcounter_tn0][pcounter_tc0 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}
										compute();
										cycle++;
										//printf(" --------------------------------------------current cycle is---%10d-----------------------------------------------\n", cycle);
										valid_vec2 = true;
										prev_pcounter_tc0 = pcounter_tc0;
										prev_pcounter_tk0 = pcounter_tk0;
										prev_pcounter_tn0 = pcounter_tn0;
									}
								}
							}
						}
						if (valid_vec2 && prev_pcounter_tn0 != -1 && prev_pcounter_tc0 != -1)
						{
							for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - prev_pcounter_tc0); i++)
							{
								SUBMAT_BO[prev_pcounter_tn0][prev_pcounter_tc0 + i] = REG_VEC2[i];
							}
							valid_vec2 = false;
							num_sram_access += MULTIPLIER_NUM;
						}

						compute_latency = cycle - prev_cycle; //recore the latency of computation
						if (dram_access_latency > compute_latency)
							cycle = prev_cycle + dram_access_latency; //if dram latency costs more cycles, update current cycle using dram latency

						prev_pcounter_k0 = pcounter_k0;
						valid_submat_bo = true;
						prev_pcounter_c0 = pcounter_c0;
						prev_pcounter_n0 = pcounter_n0;
					}
				}
			}
			//  store SUBMAT_BO
			// if SUBMAT_BO is valid, which means it has to write back to DRAM
			if (valid_submat_bo && prev_pcounter_n0 != -1 && prev_pcounter_c0 != -1)
			{
				for (int i = 0; i < std::min(TILE_SIZE_TN0, DIM_N - prev_pcounter_n0); i++)
				{
					for (int j = 0; j < std::min(TILE_SIZE_TC0, DIM_C - prev_pcounter_c0); j++)
					{
						MAT_B[i + prev_pcounter_n0][prev_pcounter_c0 + j] = SUBMAT_BO[i][j];
					}
				}
				//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - prev_pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - prev_pcounter_c0);
				num_dram_access += countDramAccess(SUBMAT_BO);
			}
			valid_submat_bo = false;

			//off-chip communication: SPMM2, O=A*B
			//load matrix A, B, O to SUBMAT_A, SUBMAT_B, SUBMAT_O
			for (pcounter_m1 = 0; pcounter_m1 < DIM_M; pcounter_m1 += TILE_SIZE_TM1)
			{
				for (pcounter_c1 = 0; pcounter_c1 < DIM_C; pcounter_c1 += TILE_SIZE_TC1)
				{
					for (pcounter_n1 = 0; pcounter_n1 < DIM_N; pcounter_n1 += TILE_SIZE_TN1)
					{
						prev_num_dram_access = num_dram_access; //record the current number of dram accesses

						if (pcounter_m1 != prev_pcounter_m1 || pcounter_n1 != prev_pcounter_n1)
						{
							//  load SUBMAT_A
							for (int i = 0; i < std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1); i++)
							{
								SUBMAT_A[i].assign(MAT_A[i + pcounter_m1].begin() + pcounter_n1, MAT_A[i + pcounter_m1].begin() + std::min(pcounter_n1 + TILE_SIZE_TN1, DIM_N));
								//SUBMAT_A[i].resize(TILE_SIZE_TN1);
							}
							//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TN1, DIM_N - pcounter_n1);
							num_dram_access += countDramAccess(SUBMAT_A);
						}

						if (pcounter_n1 != prev_pcounter_n1 || pcounter_c1 != prev_pcounter_c1)
						{
							//  load SUBMAT_BI
							for (int i = 0; i < std::min(TILE_SIZE_TN1, DIM_N - pcounter_n1); i++)
							{
								SUBMAT_BI[i].assign(MAT_B[i + pcounter_n1].begin() + pcounter_c1, MAT_B[i + pcounter_n1].begin() + std::min(pcounter_c1 + TILE_SIZE_TC1, DIM_C));
								//SUBMAT_BI[i].resize(TILE_SIZE_TC1);
							}
							//num_dram_access += std::min(TILE_SIZE_TN1, DIM_N - pcounter_n1) * std::min(TILE_SIZE_TC1, DIM_C - pcounter_c1);
							num_dram_access += countDramAccess(SUBMAT_BI);
						}

						if (pcounter_m1 != prev_pcounter_m1 || pcounter_c1 != prev_pcounter_c1)
						{
							//  store SUBMAT_O
							// if SUBMAT_BO is valid, which means it has to write back to DRAM
							if (valid_submat_o && prev_pcounter_m1 != -1 && prev_pcounter_c1 != -1)
							{
								for (int i = 0; i < std::min(TILE_SIZE_TM1, DIM_M - prev_pcounter_m1); i++)
								{
									for (int j = 0; j < std::min(TILE_SIZE_TC1, DIM_C - prev_pcounter_c1); j++)
									{
										MAT_O[i + prev_pcounter_m1][j + prev_pcounter_c1] = SUBMAT_O[i][j];
									}
								}
								//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TC1, DIM_C - pcounter_c1);
								num_dram_access += countDramAccess(SUBMAT_O);
								valid_submat_o = false;
							}

							//  load SUBMAT_O
							for (int i = 0; i < std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1); i++)
							{
								SUBMAT_O[i].assign(MAT_O[i + pcounter_m1].begin() + pcounter_c1, MAT_O[i + pcounter_m1].begin() + std::min(pcounter_c1 + TILE_SIZE_TC1, DIM_C));
								//SUBMAT_O[i].resize(TILE_SIZE_TC1);
							}
							//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TC1, DIM_C - pcounter_c1);
							num_dram_access += countDramAccess(SUBMAT_O);
						}

						dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access - prev_num_dram_access) / (DRAM_BANDWIDTH / CLOCK_RATE); //calculate the dram access latency

						//on-chip computation
						// loop order tn0->tk->tc0
						prev_cycle = cycle;
						for (pcounter_tm1 = 0; pcounter_tm1 < std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1); pcounter_tm1++)
						{
							for (pcounter_tn1 = 0; pcounter_tn1 < std::min(TILE_SIZE_TN1, DIM_N - pcounter_n1); pcounter_tn1++)
							{
								if (SUBMAT_A[pcounter_tm1][pcounter_tn1] != 0)
								{
									if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tn1 != prev_pcounter_tn1)
									{
										REG_SCALAR = SUBMAT_A[pcounter_tm1][pcounter_tn1];
										num_sram_access += 1;
									}

									for (pcounter_tc1 = 0; pcounter_tc1 < std::min(TILE_SIZE_TC1, DIM_C - pcounter_c1); pcounter_tc1 += MULTIPLIER_NUM)
									{
										// load element of BI to registers
										if (pcounter_tn1 != prev_pcounter_tn1 || pcounter_tc1 != prev_pcounter_tc1)
										{
											for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC1 - pcounter_tc1); i++)
											{
												REG_VEC1[i] = SUBMAT_BI[pcounter_tn1][pcounter_tc1 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}

										// load elements of O to registers
										if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tc1 != prev_pcounter_tc1)
										{
											// if replaced, the results in registers need write back
											if (valid_vec2 && prev_pcounter_tm1 != -1 && prev_pcounter_tc1 != -1)
											{
												for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC1 - prev_pcounter_tc1); i++)
												{
													SUBMAT_O[prev_pcounter_tm1][prev_pcounter_tc1 + i] = REG_VEC2[i];
												}
												valid_vec2 = false;
												num_sram_access += MULTIPLIER_NUM;
											}

											for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC1 - pcounter_tc1); i++)
											{
												REG_VEC2[i] = SUBMAT_O[pcounter_tm1][pcounter_tc1 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}
										compute();
										cycle++;
										//printf(" --------------------------------------------current cycle is---%10d-----------------------------------------------\n", cycle);
										valid_vec2 = true;
										prev_pcounter_tc1 = pcounter_tc1;
										prev_pcounter_tn1 = pcounter_tn1;
										prev_pcounter_tm1 = pcounter_tm1;
									}
								}
							}
						}
						if (valid_vec2 && prev_pcounter_tm1 != -1 && prev_pcounter_tc1 != -1)
						{
							for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC1 - prev_pcounter_tc1); i++)
							{
								SUBMAT_O[prev_pcounter_tm1][prev_pcounter_tc1 + i] = REG_VEC2[i];
							}
							valid_vec2 = false;
							num_sram_access += MULTIPLIER_NUM;
						}

						compute_latency = cycle - prev_cycle; //recore the latency of computation
						if (dram_access_latency > compute_latency)
							cycle = prev_cycle + dram_access_latency; //if dram latency costs more cycles, update current cycle using dram latency

						valid_submat_o = true;
						prev_pcounter_n1 = pcounter_n1;
						prev_pcounter_c1 = pcounter_c1;
						prev_pcounter_m1 = pcounter_m1;
					}
				}
			}

			//  store SUBMAT_O
			// if SUBMAT_O is valid, which means it has to write back to DRAM
			if (valid_submat_o && prev_pcounter_m1 != -1 && prev_pcounter_c1 != -1)
			{
				for (int i = 0; i < std::min(TILE_SIZE_TM1, DIM_M - prev_pcounter_m1); i++)
				{
					for (int j = 0; j < std::min(TILE_SIZE_TC1, DIM_C - prev_pcounter_c1); j++)
					{
						MAT_O[i + prev_pcounter_m1][prev_pcounter_c1 + j] = SUBMAT_O[i][j];
					}
				}
				//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - prev_pcounter_m) * std::min(TILE_SIZE_TC1, DIM_C - prev_pcounter_c1);
				num_dram_access += countDramAccess(SUBMAT_O);
			}
			valid_submat_o = false;
		}
		// Enable loop fusion
		else
		{
			prev_pcounter_n0 = -1;
			prev_pcounter_c0 = -1;
			prev_pcounter_k0 = -1;
			prev_pcounter_tn0 = -1;
			prev_pcounter_tc0 = -1;
			prev_pcounter_tk0 = -1;

			prev_pcounter_m1 = -1;
			prev_pcounter_c1 = -1;
			prev_pcounter_n1 = -1;
			prev_pcounter_tm1 = -1;
			prev_pcounter_tc1 = -1;
			prev_pcounter_tn1 = -1;

			for (pcounter_n0 = 0; pcounter_n0 < DIM_N; pcounter_n0 += TILE_SIZE_TN0)
			{
				for (pcounter_c0 = 0; pcounter_c0 < DIM_C; pcounter_c0 += TILE_SIZE_TC0)
				{
					for (pcounter_k0 = 0; pcounter_k0 < DIM_K; pcounter_k0 += TILE_SIZE_TK0)
					{
						prev_num_dram_access = num_dram_access; //record the current number of dram accesses

						if (pcounter_n0 != prev_pcounter_n0 || pcounter_k0 != prev_pcounter_k0)
						{
							//  load SUBMAT_X
							for (int i = 0; i < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); i++)
							{
								SUBMAT_X[i].assign(MAT_X[i + pcounter_n0].begin() + pcounter_k0, MAT_X[i + pcounter_n0].begin() + std::min(pcounter_k0 + TILE_SIZE_TK0, DIM_K));
								//SUBMAT_X[i].resize(TILE_SIZE_TK);
							}
							//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TK, DIM_K - pcounter_k);
							num_dram_access += countDramAccess(SUBMAT_X);
						}

						if (pcounter_k0 != prev_pcounter_k0 || pcounter_c0 != prev_pcounter_c0)
						{
							//  load SUBMAT_W
							for (int i = 0; i < std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0); i++)
							{
								SUBMAT_W[i].assign(MAT_W[i + pcounter_k0].begin() + pcounter_c0, MAT_W[i + pcounter_k0].begin() + std::min(pcounter_c0 + TILE_SIZE_TC0, DIM_C));
								//SUBMAT_W[i].resize(TILE_SIZE_TC0);
							}
							//num_dram_access += std::min(TILE_SIZE_TK, DIM_K - pcounter_k) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							num_dram_access += countDramAccess(SUBMAT_W);
						}

						if (pcounter_n0 != prev_pcounter_n0 || pcounter_c0 != prev_pcounter_c0)
						{
							//  load SUBMAT_BO
							// if SUBMAT_BO is valid, which means it has to write back to DRAM
							if (valid_submat_bo && prev_pcounter_n0 != -1 && prev_pcounter_c0 != -1)
							{
								for (int i = 0; i < std::min(TILE_SIZE_TN0, DIM_N - prev_pcounter_n0); i++)
								{
									for (int j = 0; j < std::min(TILE_SIZE_TC0, DIM_C - prev_pcounter_c0); j++)
									{
										MAT_B[i + prev_pcounter_n0][j + prev_pcounter_c0] = SUBMAT_BO[i][j];
									}
								}
								//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
								//num_dram_access += countDramAccess(SUBMAT_X);
								valid_submat_bo = false;
							}

							//  load SUBMAT_BO
							for (int i = 0; i < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); i++)
							{
								SUBMAT_BO[i].assign(MAT_B[i + pcounter_n0].begin() + pcounter_c0, MAT_B[i + pcounter_n0].begin() + std::min(pcounter_c0 + TILE_SIZE_TC0, DIM_C));
								//SUBMAT_BO[i].resize(TILE_SIZE_TC0);
							}
							//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
						}

						dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access - prev_num_dram_access) / (DRAM_BANDWIDTH / CLOCK_RATE); //calculate the dram access latency

						//on-chip computation
						// loop order tn0->tk->tc0
						for (pcounter_tn0 = 0; pcounter_tn0 < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); pcounter_tn0++)
						{
							for (pcounter_tk0 = 0; pcounter_tk0 < std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0); pcounter_tk0++)
							{
								if (SUBMAT_X[pcounter_tn0][pcounter_tk0] != 0)
								{
									if (pcounter_tn0 != prev_pcounter_tn0 || pcounter_tk0 != prev_pcounter_tk0)
									{
										REG_SCALAR = SUBMAT_X[pcounter_tn0][pcounter_tk0];
										num_sram_access += 1;
									}

									for (pcounter_tc0 = 0; pcounter_tc0 < std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0); pcounter_tc0 += MULTIPLIER_NUM)
									{
										// load element of W to registers
										if (pcounter_tk0 != prev_pcounter_tk0 || pcounter_tc0 != prev_pcounter_tc0)
										{
											for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc0); i++)
											{
												REG_VEC1[i] = SUBMAT_W[pcounter_tk0][pcounter_tc0 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}

										// load elements of B to registers
										if (pcounter_tn0 != prev_pcounter_tn0 || pcounter_tc0 != prev_pcounter_tc0)
										{
											// if replaced, the results in registers need write back
											if (valid_vec2 && prev_pcounter_tn0 != -1 && prev_pcounter_tc0 != -1)
											{
												for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - prev_pcounter_tc0); i++)
												{
													SUBMAT_BO[prev_pcounter_tn0][prev_pcounter_tc0 + i] = REG_VEC2[i];
												}
												valid_vec2 = false;
												num_sram_access += MULTIPLIER_NUM;
											}

											for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc0); i++)
											{
												REG_VEC2[i] = SUBMAT_BO[pcounter_tn0][pcounter_tc0 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}
										compute();
										cycle++;

										//printf(" --------------------------------------------current cycle is---%10d-----------------------------------------------\n", cycle);
										valid_vec2 = true;
										prev_pcounter_tc0 = pcounter_tc0;
										prev_pcounter_tk0 = pcounter_tk0;
										prev_pcounter_tn0 = pcounter_tn0;
									}
								}
							}
						}
						if (valid_vec2 && prev_pcounter_tn0 != -1 && prev_pcounter_tc0 != -1)
						{
							for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - prev_pcounter_tc0); i++)
							{
								SUBMAT_BO[prev_pcounter_tn0][prev_pcounter_tc0 + i] = REG_VEC2[i];
							}
							valid_vec2 = false;
							num_sram_access += MULTIPLIER_NUM;
						}

						compute_latency = cycle - prev_cycle; //recore the latency of computation
						if (dram_access_latency > compute_latency)
							cycle = prev_cycle + dram_access_latency; //if dram latency costs more cycles, update current cycle using dram latency

						prev_pcounter_k0 = pcounter_k0;
						valid_submat_bo = true;
						prev_pcounter_c0 = pcounter_c0;
						prev_pcounter_n0 = pcounter_n0;
					}

					for (pcounter_m1 = 0; pcounter_m1 < DIM_M; pcounter_m1 += TILE_SIZE_TM1)
					{
						prev_num_dram_access = num_dram_access; //record the current number of dram accesses

						// off-chip communication, load SUBMAT_A, SUBMAT_O
						if (pcounter_m1 != prev_pcounter_m1 || pcounter_n0 != prev_pcounter_n0)
						{
							//  load SUBMAT_A
							for (int i = 0; i < std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1); i++)
							{
								SUBMAT_A[i].assign(MAT_A[i + pcounter_m1].begin() + pcounter_n0, MAT_A[i + pcounter_m1].begin() + std::min(pcounter_n0 + TILE_SIZE_TN0, DIM_N));
								//SUBMAT_A[i].resize(TILE_SIZE_TN1);
							}
							//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
							num_dram_access += countDramAccess(SUBMAT_A);
						}

						if (pcounter_n0 != prev_pcounter_n0 || pcounter_c0 != prev_pcounter_c0)
						{
							//  load SUBMAT_BI
							for (int i = 0; i < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); i++)
							{
								SUBMAT_BI[i].assign(MAT_B[i + pcounter_n0].begin() + pcounter_c0, MAT_B[i + pcounter_n0].begin() + std::min(pcounter_c0 + TILE_SIZE_TC0, DIM_C));
								//SUBMAT_BI[i].resize(TILE_SIZE_TC1);
							}
							//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							num_dram_access += countDramAccess(SUBMAT_BI);
						}

						if (pcounter_m1 != prev_pcounter_m1 || pcounter_c0 != prev_pcounter_c0)
						{
							//  store SUBMAT_O
							// if SUBMAT_BO is valid, which means it has to write back to DRAM
							if (valid_submat_o && prev_pcounter_m1 != -1 && prev_pcounter_c0 != -1)
							{
								for (int i = 0; i < std::min(TILE_SIZE_TM1, DIM_M - prev_pcounter_m1); i++)
								{
									for (int j = 0; j < std::min(TILE_SIZE_TC0, DIM_C - prev_pcounter_c0); j++)
									{
										MAT_O[i + prev_pcounter_m1][j + prev_pcounter_c0] = SUBMAT_O[i][j];
									}
								}
								//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
								num_dram_access += countDramAccess(SUBMAT_O);
								valid_submat_o = false;
							}

							//  load SUBMAT_O
							for (int i = 0; i < std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1); i++)
							{
								SUBMAT_O[i].assign(MAT_O[i + pcounter_m1].begin() + pcounter_c0, MAT_O[i + pcounter_m1].begin() + std::min(pcounter_c0 + TILE_SIZE_TC0, DIM_C));
								//SUBMAT_O[i].resize(TILE_SIZE_TC1);
							}
							//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							num_dram_access += countDramAccess(SUBMAT_O);
						}

						dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access - prev_num_dram_access) / (DRAM_BANDWIDTH / CLOCK_RATE); //calculate the dram access latency

						//on-chip computation
						// loop order tn0->tk->tc0
						for (pcounter_tm1 = 0; pcounter_tm1 < std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1); pcounter_tm1++)
						{
							for (pcounter_tn1 = 0; pcounter_tn1 < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); pcounter_tn1++)
							{
								if (SUBMAT_A[pcounter_tm1][pcounter_tn1] != 0)
								{
									if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tn1 != prev_pcounter_tn1)
									{
										REG_SCALAR = SUBMAT_A[pcounter_tm1][pcounter_tn1];
										num_sram_access += 1;
									}

									for (pcounter_tc1 = 0; pcounter_tc1 < std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0); pcounter_tc1 += MULTIPLIER_NUM)
									{
										// load element of BI to registers
										if (pcounter_tn1 != prev_pcounter_tn1 || pcounter_tc1 != prev_pcounter_tc1)
										{
											for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc1); i++)
											{
												REG_VEC1[i] = SUBMAT_BO[pcounter_tn1][pcounter_tc1 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}

										// load elements of O to registers
										if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tc1 != prev_pcounter_tc1)
										{
											// if replaced, the results in registers need write back
											if (valid_vec2 && prev_pcounter_tm1 != -1 && prev_pcounter_tc1 != -1)
											{
												for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - prev_pcounter_tc1); i++)
												{
													SUBMAT_O[prev_pcounter_tm1][prev_pcounter_tc1 + i] = REG_VEC2[i];
												}
												valid_vec2 = false;
												num_sram_access += MULTIPLIER_NUM;
											}

											for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc1); i++)
											{
												REG_VEC2[i] = SUBMAT_O[pcounter_tm1][pcounter_tc1 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}
										compute();
										cycle++;
										//printf(" --------------------------------------------current cycle is---%10d-----------------------------------------------\n", cycle);
										valid_vec2 = true;
										prev_pcounter_tc1 = pcounter_tc1;
										prev_pcounter_tn1 = pcounter_tn1;
										prev_pcounter_tm1 = pcounter_tm1;
									}
								}
							}
						}
						if (valid_vec2 && prev_pcounter_tm1 != -1 && prev_pcounter_tc1 != -1)
						{
							for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - prev_pcounter_tc1); i++)
							{
								SUBMAT_O[prev_pcounter_tm1][prev_pcounter_tc1 + i] = REG_VEC2[i];
							}
							valid_vec2 = false;
							num_sram_access += MULTIPLIER_NUM;
						}

						compute_latency = cycle - prev_cycle; //recore the latency of computation
						if (dram_access_latency > compute_latency)
							cycle = prev_cycle + dram_access_latency; //if dram latency costs more cycles, update current cycle using dram latency

						valid_submat_o = true;
						prev_pcounter_n1 = pcounter_n1;
						prev_pcounter_c1 = pcounter_c1;
						prev_pcounter_m1 = pcounter_m1;
						prev_pcounter_n0 = pcounter_n0;
						prev_pcounter_c0 = pcounter_c0;
						prev_pcounter_m0 = pcounter_m0;
					}
				}
			}
			//  store SUBMAT_O
			// if SUBMAT_O is valid, which means it has to write back to DRAM
			if (valid_submat_o && prev_pcounter_m1 != -1 && prev_pcounter_c0 != -1)
			{
				for (int i = 0; i < std::min(TILE_SIZE_TM1, DIM_M - prev_pcounter_m1); i++)
				{
					for (int j = 0; j < std::min(TILE_SIZE_TC0, DIM_C - prev_pcounter_c0); j++)
					{
						MAT_O[i + prev_pcounter_m1][prev_pcounter_c0 + j] = SUBMAT_O[i][j];
					}
				}
				//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - prev_pcounter_m) * std::min(TILE_SIZE_TC0, DIM_C - prev_pcounter_c1);
				num_dram_access += countDramAccess(SUBMAT_O);
			}
			valid_submat_o = false;
		}
	}
	else if (ACCELERATOR == "HYGCN")
	{
		// Notes: the exeuction order of HyGCN is different from that of GCNAX/AWBGCN
		// HyGCN, B[M,K]=A[M,N]*X[N,K], O[M,C]=B[M, K]*W[K, C]; The dimensions of B is different,

		prev_pcounter_m0 = -1;
		prev_pcounter_k0 = -1;
		prev_pcounter_n0 = -1;

		prev_pcounter_tm0 = -1;
		prev_pcounter_tk0 = -1;
		prev_pcounter_tn0 = -1;

		prev_pcounter_m1 = -1;
		prev_pcounter_k1 = -1;
		prev_pcounter_c1 = -1;

		prev_pcounter_tm1 = -1;
		prev_pcounter_tk1 = -1;
		prev_pcounter_tc1 = -1;

		for (pcounter_m0 = 0; pcounter_m0 < DIM_M; pcounter_m0 += TILE_SIZE_TM0)
		{
			for (pcounter_k0 = 0; pcounter_k0 < DIM_K; pcounter_k0 += TILE_SIZE_TK0)
			{
				for (pcounter_n0 = 0; pcounter_n0 < DIM_N; pcounter_n0 += TILE_SIZE_TN0)
				{
					prev_num_dram_access = num_dram_access; //record the current number of dram accesses
					if (pcounter_m0 != prev_pcounter_m0 || pcounter_n0 != prev_pcounter_n0)
					{
						//load SUBMAT_A
						for (int i = 0; i < std::min(TILE_SIZE_TM0, DIM_M - pcounter_m0); i++)
						{
							SUBMAT_A[i].assign(MAT_A[i + pcounter_m0].begin() + pcounter_n0, MAT_A[i + pcounter_m0].begin() + std::min(pcounter_n0 + TILE_SIZE_TN0, DIM_N));
							//SUBMAT_X[i].resize(TILE_SIZE_TK);
						}
						//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TK, DIM_K - pcounter_k);
						num_dram_access += countDramAccess(SUBMAT_A);
					}

					if (pcounter_n0 != prev_pcounter_n0 || pcounter_k0 != prev_pcounter_k0)
					{
						//  load SUBMAT_W
						for (int i = 0; i < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); i++)
						{
							SUBMAT_X[i].assign(MAT_X[i + pcounter_n0].begin() + pcounter_k0, MAT_X[i + pcounter_n0].begin() + std::min(pcounter_k0 + TILE_SIZE_TK0, DIM_K));
							//SUBMAT_W[i].resize(TILE_SIZE_TC0);
						}
						//num_dram_access += std::min(TILE_SIZE_TK, DIM_K - pcounter_k) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
						num_dram_access += countDramAccess(SUBMAT_X);
					}

					if (pcounter_m0 != prev_pcounter_m0 || pcounter_k0 != prev_pcounter_k0)
					{
						//  load SUBMAT_BO
						// if SUBMAT_BO is valid, which means it has to write back to DRAM
						if (valid_submat_bo && prev_pcounter_m0 != -1 && prev_pcounter_k0 != -1)
						{
							for (int i = 0; i < std::min(TILE_SIZE_TM0, DIM_M - prev_pcounter_m0); i++)
							{
								for (int j = 0; j < std::min(TILE_SIZE_TK0, DIM_K - prev_pcounter_k0); j++)
								{
									MAT_B[i + prev_pcounter_m0][j + prev_pcounter_k0] = SUBMAT_BO[i][j];
								}
							}
							//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							//num_dram_access += countDramAccess(SUBMAT_X);
							valid_submat_bo = false;
						}

						//  load SUBMAT_BO
						for (int i = 0; i < std::min(TILE_SIZE_TM0, DIM_M - pcounter_m0); i++)
						{
							SUBMAT_BO[i].assign(MAT_B[i + pcounter_m0].begin() + pcounter_k0, MAT_B[i + pcounter_m0].begin() + std::min(pcounter_k0 + TILE_SIZE_TK0, DIM_K));
							//SUBMAT_BO[i].resize(TILE_SIZE_TC0);
						}
						//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
					}

					dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access - prev_num_dram_access) / (DRAM_BANDWIDTH / CLOCK_RATE); //calculate the dram access latency

					//on-chip computation
					// loop order tn0->tk->tc0
					for (pcounter_tm0 = 0; pcounter_tm0 < std::min(TILE_SIZE_TM0, DIM_M - pcounter_m0); pcounter_tm0++)
					{
						for (pcounter_tn0 = 0; pcounter_tn0 < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); pcounter_tn0++)
						{
							if (SUBMAT_A[pcounter_tm0][pcounter_tn0] != 0)
							{
								if (pcounter_tm0 != prev_pcounter_tm0 || pcounter_tn0 != prev_pcounter_tn0)
								{
									REG_SCALAR = SUBMAT_A[pcounter_tm0][pcounter_tn0];
									num_sram_access += 1;
								}

								for (pcounter_tk0 = 0; pcounter_tk0 < std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0); pcounter_tk0 += MULTIPLIER_NUM)
								{
									// load element of W to registers
									if (pcounter_tn0 != prev_pcounter_tn0 || pcounter_tk0 != prev_pcounter_tk0)
									{
										for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TK0 - pcounter_tk0); i++)
										{
											REG_VEC1[i] = SUBMAT_X[pcounter_tn0][pcounter_tk0 + i];
										}
										num_sram_access += MULTIPLIER_NUM;
									}

									// load elements of B to registers
									if (pcounter_tm0 != prev_pcounter_tm0 || pcounter_tk0 != prev_pcounter_tk0)
									{
										// if replaced, the results in registers need write back
										if (valid_vec2 && prev_pcounter_tm0 != -1 && prev_pcounter_tk0 != -1)
										{
											for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TK0 - prev_pcounter_tk0); i++)
											{
												SUBMAT_BO[prev_pcounter_tm0][prev_pcounter_tk0 + i] = REG_VEC2[i];
											}
											valid_vec2 = false;
											num_sram_access += MULTIPLIER_NUM;
										}

										for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TK0 - pcounter_tk0); i++)
										{
											REG_VEC2[i] = SUBMAT_BO[pcounter_tm0][pcounter_tk0 + i];
										}
										num_sram_access += MULTIPLIER_NUM;
									}
									compute();
									cycle++;
									//printf(" --------------------------------------------current cycle is---%10d-----------------------------------------------\n", cycle);
									valid_vec2 = true;
									prev_pcounter_tk0 = pcounter_tk0;
									prev_pcounter_tn0 = pcounter_tn0;
									prev_pcounter_tm0 = pcounter_tm0;
								}
							}
						}
					}
					if (valid_vec2 && prev_pcounter_tm0 != -1 && prev_pcounter_tk0 != -1)
					{
						for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TK0 - prev_pcounter_tk0); i++)
						{
							SUBMAT_BO[prev_pcounter_tm0][prev_pcounter_tk0 + i] = REG_VEC2[i];
						}
						valid_vec2 = false;
						num_sram_access += MULTIPLIER_NUM;
					}

					compute_latency = cycle - prev_cycle; //recore the latency of computation
					if (dram_access_latency > compute_latency)
						cycle = prev_cycle + dram_access_latency; //if dram latency costs more cycles, update current cycle using dram latency

					prev_pcounter_n0 = pcounter_n0;
					valid_submat_bo = true;
					prev_pcounter_k0 = pcounter_k0;
					prev_pcounter_m0 = pcounter_m0;
				}

				for (pcounter_c1 = 0; pcounter_c1 < DIM_C; pcounter_c1 += TILE_SIZE_TC1)
				{
					prev_num_dram_access = num_dram_access; //record the current number of dram accesses

					// off-chip communication, load SUBMAT_A, SUBMAT_O
					if (pcounter_c1 != prev_pcounter_c1 || pcounter_k0 != prev_pcounter_k0)
					{
						//  load SUBMAT_A
						for (int i = 0; i < std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0); i++)
						{
							SUBMAT_W[i].assign(MAT_W[i + pcounter_k0].begin() + pcounter_c1, MAT_W[i + pcounter_k0].begin() + std::min(pcounter_c1 + TILE_SIZE_TC1, DIM_C));
							//SUBMAT_A[i].resize(TILE_SIZE_TN1);
						}
						//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
						num_dram_access += countDramAccess(SUBMAT_W);
					}

					if (pcounter_m0 != prev_pcounter_m0 || pcounter_k0 != prev_pcounter_k0)
					{
						//  load SUBMAT_BI
						for (int i = 0; i < std::min(TILE_SIZE_TM0, DIM_M - pcounter_m0); i++)
						{
							SUBMAT_BI[i].assign(MAT_B[i + pcounter_m0].begin() + pcounter_k0, MAT_B[i + pcounter_m0].begin() + std::min(pcounter_k0 + TILE_SIZE_TK0, DIM_K));
							//SUBMAT_BI[i].resize(TILE_SIZE_TC1);
						}
						//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
						num_dram_access += countDramAccess(SUBMAT_BI);
					}

					if (pcounter_m0 != prev_pcounter_m0 || pcounter_c1 != prev_pcounter_c1)
					{
						//  store SUBMAT_O
						// if SUBMAT_BO is valid, which means it has to write back to DRAM
						if (valid_submat_o && prev_pcounter_m0 != -1 && prev_pcounter_c1 != -1)
						{
							for (int i = 0; i < std::min(TILE_SIZE_TM0, DIM_M - prev_pcounter_m0); i++)
							{
								for (int j = 0; j < std::min(TILE_SIZE_TC1, DIM_C - prev_pcounter_c1); j++)
								{
									MAT_O[i + prev_pcounter_m0][j + prev_pcounter_c1] = SUBMAT_O[i][j];
								}
							}
							//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							num_dram_access += countDramAccess(SUBMAT_O);
							valid_submat_o = false;
						}

						//  load SUBMAT_O
						for (int i = 0; i < std::min(TILE_SIZE_TM0, DIM_M - pcounter_m0); i++)
						{
							SUBMAT_O[i].assign(MAT_O[i + pcounter_m0].begin() + pcounter_c1, MAT_O[i + pcounter_m0].begin() + std::min(pcounter_c1 + TILE_SIZE_TC1, DIM_C));
							//SUBMAT_O[i].resize(TILE_SIZE_TC1);
						}
						//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
						num_dram_access += countDramAccess(SUBMAT_O);
					}

					dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access - prev_num_dram_access) / (DRAM_BANDWIDTH / CLOCK_RATE); //calculate the dram access latency

					//on-chip computation
					// loop order tn0->tk->tc0
					for (pcounter_tm1 = 0; pcounter_tm1 < std::min(TILE_SIZE_TM0, DIM_M - pcounter_m0); pcounter_tm1++)
					{
						for (pcounter_tk1 = 0; pcounter_tk1 < std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0); pcounter_tk1++)
						{
							if (SUBMAT_BO[pcounter_tm1][pcounter_tk1] != 0)
							{
								if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tk1 != prev_pcounter_tk1)
								{
									REG_SCALAR = SUBMAT_BO[pcounter_tm1][pcounter_tk1];
									num_sram_access += 1;
								}

								for (pcounter_tc1 = 0; pcounter_tc1 < std::min(TILE_SIZE_TC1, DIM_C - pcounter_c1); pcounter_tc1 += MULTIPLIER_NUM)
								{
									// load element of BI to registers
									if (pcounter_tk1 != prev_pcounter_tk1 || pcounter_tc1 != prev_pcounter_tc1)
									{
										for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC1 - pcounter_tc1); i++)
										{
											REG_VEC1[i] = SUBMAT_W[pcounter_tk1][pcounter_tc1 + i];
										}
										num_sram_access += MULTIPLIER_NUM;
									}

									// load elements of O to registers
									if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tc1 != prev_pcounter_tc1)
									{
										// if replaced, the results in registers need write back
										if (valid_vec2 && prev_pcounter_tm1 != -1 && prev_pcounter_tc1 != -1)
										{
											for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC1 - prev_pcounter_tc1); i++)
											{
												SUBMAT_O[prev_pcounter_tm1][prev_pcounter_tc1 + i] = REG_VEC2[i];
											}
											valid_vec2 = false;
											num_sram_access += MULTIPLIER_NUM;
										}

										for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC1 - pcounter_tc1); i++)
										{
											REG_VEC2[i] = SUBMAT_O[pcounter_tm1][pcounter_tc1 + i];
										}
										num_sram_access += MULTIPLIER_NUM;
									}
									compute();
									cycle++;
									//printf(" --------------------------------------------current cycle is---%10d-----------------------------------------------\n", cycle);
									valid_vec2 = true;
									prev_pcounter_tc1 = pcounter_tc1;
									prev_pcounter_tk1 = pcounter_tk1;
									prev_pcounter_tm1 = pcounter_tm1;
								}
							}
						}
					}
					if (valid_vec2 && prev_pcounter_tm1 != -1 && prev_pcounter_tc1 != -1)
					{
						for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC1 - prev_pcounter_tc1); i++)
						{
							SUBMAT_O[prev_pcounter_tm1][prev_pcounter_tc1 + i] = REG_VEC2[i];
						}
						valid_vec2 = false;
						num_sram_access += MULTIPLIER_NUM;
					}

					compute_latency = cycle - prev_cycle; //recore the latency of computation
					if (dram_access_latency > compute_latency)
						cycle = prev_cycle + dram_access_latency; //if dram latency costs more cycles, update current cycle using dram latency

					valid_submat_o = true;
					prev_pcounter_k1 = pcounter_k1;
					prev_pcounter_c1 = pcounter_c1;
					prev_pcounter_m1 = pcounter_m1;
					prev_pcounter_k0 = pcounter_k0;
					prev_pcounter_c0 = pcounter_c0;
					prev_pcounter_m0 = pcounter_m0;
				}
			}
		}
		//  store SUBMAT_O
		// if SUBMAT_O is valid, which means it has to write back to DRAM
		if (valid_submat_o && prev_pcounter_m0 != -1 && prev_pcounter_c1 != -1)
		{
			for (int i = 0; i < std::min(TILE_SIZE_TM0, DIM_M - prev_pcounter_m0); i++)
			{
				for (int j = 0; j < std::min(TILE_SIZE_TC1, DIM_C - prev_pcounter_c1); j++)
				{
					MAT_O[i + prev_pcounter_m0][prev_pcounter_c1 + j] = SUBMAT_O[i][j];
				}
			}
			//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - prev_pcounter_m) * std::min(TILE_SIZE_TC0, DIM_C - prev_pcounter_c1);
			num_dram_access += countDramAccess(SUBMAT_O);
		}
		valid_submat_o = false;
	}
	else if (ACCELERATOR == "AWBGCN")
	{
		COMPUTE_TYPE = 0;

		prev_pcounter_n0 = -1;
		prev_pcounter_c0 = -1;
		prev_pcounter_k0 = -1;
		prev_pcounter_tn0 = -1;
		prev_pcounter_tc0 = -1;
		prev_pcounter_tk0 = -1;

		prev_pcounter_m1 = -1;
		prev_pcounter_c1 = -1;
		prev_pcounter_n1 = -1;
		prev_pcounter_tm1 = -1;
		prev_pcounter_tc1 = -1;
		prev_pcounter_tn1 = -1;

		for (pcounter_n0 = 0; pcounter_n0 < DIM_N; pcounter_n0 += TILE_SIZE_TN0)
		{
			for (pcounter_c0 = 0; pcounter_c0 < DIM_C; pcounter_c0 += TILE_SIZE_TC0)
			{
				for (pcounter_k0 = 0; pcounter_k0 < DIM_K; pcounter_k0 += TILE_SIZE_TK0)
				{
					prev_num_dram_access = num_dram_access; //record the current number of dram accesses
					if (pcounter_n0 != prev_pcounter_n0 || pcounter_k0 != prev_pcounter_k0)
					{
						//  load SUBMAT_X
						for (int i = 0; i < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); i++)
						{
							SUBMAT_X[i].assign(MAT_X[i + pcounter_n0].begin() + pcounter_k0, MAT_X[i + pcounter_n0].begin() + std::min(pcounter_k0 + TILE_SIZE_TK0, DIM_K));
							//SUBMAT_X[i].resize(TILE_SIZE_TK);
						}
						//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TK, DIM_K - pcounter_k);
						num_dram_access += countDramAccess(SUBMAT_X);
					}

					if (pcounter_k0 != prev_pcounter_k0 || pcounter_c0 != prev_pcounter_c0)
					{
						//  load SUBMAT_W
						for (int i = 0; i < std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0); i++)
						{
							SUBMAT_W[i].assign(MAT_W[i + pcounter_k0].begin() + pcounter_c0, MAT_W[i + pcounter_k0].begin() + std::min(pcounter_c0 + TILE_SIZE_TC0, DIM_C));
							//SUBMAT_W[i].resize(TILE_SIZE_TC0);
						}
						//num_dram_access += std::min(TILE_SIZE_TK, DIM_K - pcounter_k) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
						num_dram_access += countDramAccess(SUBMAT_W);
					}

					if (pcounter_n0 != prev_pcounter_n0 || pcounter_c0 != prev_pcounter_c0)
					{
						//  load SUBMAT_BO
						// if SUBMAT_BO is valid, which means it has to write back to DRAM
						if (valid_submat_bo && prev_pcounter_n0 != -1 && prev_pcounter_c0 != -1)
						{
							for (int i = 0; i < std::min(TILE_SIZE_TN0, DIM_N - prev_pcounter_n0); i++)
							{
								for (int j = 0; j < std::min(TILE_SIZE_TC0, DIM_C - prev_pcounter_c0); j++)
								{
									MAT_B[i + prev_pcounter_n0][j + prev_pcounter_c0] = SUBMAT_BO[i][j];
								}
							}
							//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							//num_dram_access += countDramAccess(SUBMAT_X);
							valid_submat_bo = false;
						}

						//  load SUBMAT_BO
						for (int i = 0; i < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); i++)
						{
							SUBMAT_BO[i].assign(MAT_B[i + pcounter_n0].begin() + pcounter_c0, MAT_B[i + pcounter_n0].begin() + std::min(pcounter_c0 + TILE_SIZE_TC0, DIM_C));
							//SUBMAT_BO[i].resize(TILE_SIZE_TC0);
						}
						//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
					}

					dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access - prev_num_dram_access) / (DRAM_BANDWIDTH / CLOCK_RATE); //calculate the dram access latency

					//on-chip computation
					// loop order tn0->tk->tc0
					for (pcounter_tn0 = 0; pcounter_tn0 < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); pcounter_tn0++)
					{
						for (pcounter_tc0 = 0; pcounter_tc0 < std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0); pcounter_tc0++)
						{
							if (true) // use inner product, cannot skip zero computations
							{

								if (pcounter_tn0 != prev_pcounter_tn0 || pcounter_tc0 != prev_pcounter_tc0)
								{
									// if the element has to be written back to SRAM
									if (valid_scalar && prev_pcounter_tn0 != -1 && prev_pcounter_tc0 != -1)
									{
										SUBMAT_BO[prev_pcounter_tn0][prev_pcounter_tc0] = REG_SCALAR;
									}

									REG_SCALAR = SUBMAT_BO[pcounter_tn0][pcounter_tc0];
									num_sram_access += 1;
								}

								for (pcounter_tk0 = 0; pcounter_tk0 < std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0); pcounter_tk0 += MULTIPLIER_NUM)
								{
									// load element of W to registers
									if (pcounter_tk0 != prev_pcounter_tk0 || pcounter_tc0 != prev_pcounter_tc0)
									{
										for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TK0 - pcounter_tk0); i++)
										{
											REG_VEC1[i] = SUBMAT_W[pcounter_tk0 + i][pcounter_tc0];
										}
										num_sram_access += MULTIPLIER_NUM;
									}

									// load elements of B to registers
									if (pcounter_tn0 != prev_pcounter_tn0 || pcounter_tk0 != prev_pcounter_tk0)
									{

										for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TK0 - pcounter_tk0); i++)
										{
											REG_VEC2[i] = SUBMAT_X[pcounter_tn0][pcounter_tk0 + i];
										}
										num_sram_access += MULTIPLIER_NUM;
									}

									int num_nonzero_computation = 0;
									for (int i = 0; i < MULTIPLIER_NUM; i++)
									{
										if (REG_VEC1[i] != 0 || REG_VEC2[i] != 0)
											num_nonzero_computation++;
									}

									compute();
									cycle += num_nonzero_computation;
									//printf(" --------------------------------------------current cycle is---%10d-----------------------------------------------\n", cycle);
									valid_vec2 = true;
									prev_pcounter_tc0 = pcounter_tc0;
									prev_pcounter_tk0 = pcounter_tk0;
									prev_pcounter_tn0 = pcounter_tn0;
								}
							}
						}
					}
					if (valid_scalar && prev_pcounter_tn0 != -1 && prev_pcounter_tc0 != -1)
					{
						SUBMAT_BO[prev_pcounter_tn0][prev_pcounter_tc0] = REG_SCALAR;
						valid_scalar = false;
						num_sram_access += 1;
					}

					compute_latency = cycle - prev_cycle; //recore the latency of computation
					if (dram_access_latency > compute_latency / MULTIPLIER_NUM)
						cycle = prev_cycle + dram_access_latency * MULTIPLIER_NUM; //if dram latency costs more cycles, update current cycle using dram latency

					prev_pcounter_k0 = pcounter_k0;
					valid_submat_bo = true;
					prev_pcounter_c0 = pcounter_c0;
					prev_pcounter_n0 = pcounter_n0;
				}

				for (pcounter_m1 = 0; pcounter_m1 < DIM_M; pcounter_m1 += TILE_SIZE_TM1)
				{
					prev_num_dram_access = num_dram_access; //record the current number of dram accesses

					// off-chip communication, load SUBMAT_A, SUBMAT_O
					if (pcounter_m1 != prev_pcounter_m1 || pcounter_n0 != prev_pcounter_n0)
					{
						//  load SUBMAT_A
						for (int i = 0; i < std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1); i++)
						{
							SUBMAT_A[i].assign(MAT_A[i + pcounter_m1].begin() + pcounter_n0, MAT_A[i + pcounter_m1].begin() + std::min(pcounter_n0 + TILE_SIZE_TN0, DIM_N));
							//SUBMAT_A[i].resize(TILE_SIZE_TN1);
						}
						//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
						num_dram_access += countDramAccess(SUBMAT_A);
					}

					if (pcounter_n0 != prev_pcounter_n0 || pcounter_c0 != prev_pcounter_c0)
					{
						//  load SUBMAT_BI
						for (int i = 0; i < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); i++)
						{
							SUBMAT_BI[i].assign(MAT_B[i + pcounter_n0].begin() + pcounter_c0, MAT_B[i + pcounter_n0].begin() + std::min(pcounter_c0 + TILE_SIZE_TC0, DIM_C));
							//SUBMAT_BI[i].resize(TILE_SIZE_TC1);
						}
						//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
						num_dram_access += countDramAccess(SUBMAT_BI);
					}

					if (pcounter_m1 != prev_pcounter_m1 || pcounter_c0 != prev_pcounter_c0)
					{
						//  store SUBMAT_O
						// if SUBMAT_BO is valid, which means it has to write back to DRAM
						if (valid_submat_o && prev_pcounter_m1 != -1 && prev_pcounter_c0 != -1)
						{
							for (int i = 0; i < std::min(TILE_SIZE_TM1, DIM_M - prev_pcounter_m1); i++)
							{
								for (int j = 0; j < std::min(TILE_SIZE_TC0, DIM_C - prev_pcounter_c0); j++)
								{
									MAT_O[i + prev_pcounter_m1][j + prev_pcounter_c0] = SUBMAT_O[i][j];
								}
							}
							//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							num_dram_access += countDramAccess(SUBMAT_O);
							valid_submat_o = false;
						}

						//  load SUBMAT_O
						for (int i = 0; i < std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1); i++)
						{
							SUBMAT_O[i].assign(MAT_O[i + pcounter_m1].begin() + pcounter_c0, MAT_O[i + pcounter_m1].begin() + std::min(pcounter_c0 + TILE_SIZE_TC0, DIM_C));
							//SUBMAT_O[i].resize(TILE_SIZE_TC1);
						}
						//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
						num_dram_access += countDramAccess(SUBMAT_O);
					}

					dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access - prev_num_dram_access) / (DRAM_BANDWIDTH / CLOCK_RATE); //calculate the dram access latency

					//on-chip computation
					// loop order tn0->tk->tc0
					for (pcounter_tm1 = 0; pcounter_tm1 < std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1); pcounter_tm1++)
					{
						for (pcounter_tc1 = 0; pcounter_tc1 < std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0); pcounter_tc1++)
						{
							if (true)
							{
								if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tc1 != prev_pcounter_tc1)
								{
									// if the element has to be written back to SRAM
									if (valid_scalar && prev_pcounter_tm1 != -1 && prev_pcounter_tc1 != -1)
									{
										SUBMAT_O[prev_pcounter_tm1][prev_pcounter_tc1] = REG_SCALAR;
									}
									REG_SCALAR = SUBMAT_O[pcounter_tm1][pcounter_tc1];
									num_sram_access += 1;
								}

								for (pcounter_tn1 = 0; pcounter_tn1 < std::min(TILE_SIZE_TN0, DIM_C - pcounter_n0); pcounter_tn1 += MULTIPLIER_NUM)
								{
									// load element of BI to registers
									if (pcounter_tn1 != prev_pcounter_tn1 || pcounter_tc1 != prev_pcounter_tc1)
									{
										for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TN0 - pcounter_tn1); i++)
										{
											REG_VEC1[i] = SUBMAT_BO[pcounter_tn1 + i][pcounter_tc1];
										}
										num_sram_access += MULTIPLIER_NUM;
									}

									// load elements of O to registers
									if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tn1 != prev_pcounter_tn1)
									{

										for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TN0 - pcounter_tn1); i++)
										{
											REG_VEC2[i] = SUBMAT_A[pcounter_tm1][pcounter_tn1 + i];
										}
										num_sram_access += MULTIPLIER_NUM;
									}
									int num_nonzero_computation = 0;
									for (int i = 0; i < MULTIPLIER_NUM; i++)
									{
										if (REG_VEC1[i] != 0 || REG_VEC2[i] != 0)
											num_nonzero_computation++;
									}

									compute();
									cycle += num_nonzero_computation;

									//printf(" --------------------------------------------current cycle is---%10d-----------------------------------------------\n", cycle);
									valid_vec2 = true;
									prev_pcounter_tc1 = pcounter_tc1;
									prev_pcounter_tn1 = pcounter_tn1;
									prev_pcounter_tm1 = pcounter_tm1;
								}
							}
						}
					}
					if (valid_scalar && prev_pcounter_tm1 != -1 && prev_pcounter_tc1 != -1)
					{
						SUBMAT_O[prev_pcounter_tm1][prev_pcounter_tc1] = REG_SCALAR;
						valid_scalar = false;
						num_sram_access += 1;
					}

					compute_latency = cycle - prev_cycle; //recore the latency of computation
					if (dram_access_latency > compute_latency / MULTIPLIER_NUM)
						cycle = prev_cycle + dram_access_latency * MULTIPLIER_NUM; //if dram latency costs more cycles, update current cycle using dram latency

					valid_submat_o = true;
					prev_pcounter_n1 = pcounter_n1;
					prev_pcounter_c1 = pcounter_c1;
					prev_pcounter_m1 = pcounter_m1;
					prev_pcounter_n0 = pcounter_n0;
					prev_pcounter_c0 = pcounter_c0;
					prev_pcounter_m0 = pcounter_m0;
				}
			}
		}
		//  store SUBMAT_O
		// if SUBMAT_O is valid, which means it has to write back to DRAM
		if (valid_submat_o && prev_pcounter_m1 != -1 && prev_pcounter_c0 != -1)
		{
			for (int i = 0; i < std::min(TILE_SIZE_TM1, DIM_M - prev_pcounter_m1); i++)
			{
				for (int j = 0; j < std::min(TILE_SIZE_TC0, DIM_C - prev_pcounter_c0); j++)
				{
					MAT_O[i + prev_pcounter_m1][prev_pcounter_c0 + j] = SUBMAT_O[i][j];
				}
			}
			//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - prev_pcounter_m) * std::min(TILE_SIZE_TC0, DIM_C - prev_pcounter_c1);
			num_dram_access += countDramAccess(SUBMAT_O);
		}
		valid_submat_o = false;
		cycle = cycle / MULTIPLIER_NUM;
	}
	else
	{
		std::cout << "Wrong Accelerator Type!" << std::endl;
	}

	end_sim = clock();
	double cost_time = (double)(end_sim - begin_sim);

	//save the results to result.txt
	//	writeResultToDram(DATASET);

	printf("\n**********************************************************************************************************************\n");
	printf(" ------------------------------------------------Sim-1 has done-------------------------------------------------------\n");
	printf(" ---------------------------------------Total cycle number is---%10d-----------------------------------------------\n", cycle - 1);
	printf(" --------------------------------------Number of DRAM accesses---%10d---------------------------------------------\n", num_dram_access);
	printf(" --------------------------------------Number of SRAM accesses---%10d---------------------------------------------\n", num_sram_access);
	//FILE *freport_test = fopen("test_report.txt","w+");
	//fprintf(freport_test, "\n\nexecution all time is : %lf seconds.\n\n", cost_time / cycle_TCK);
	//fclose(freport_test);
	printf("**********************************************************************************************************************\n");
}

void generateVirtualData(std::vector<std::vector<double>> &mat, int row, int col, double density)
{
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			double val = ((double)(std::rand() % 10000)) / 10000;
			if (val < density)
			{
				mat[i][j] = ((double)(std::rand() % 10000)) / 10000;
			}
			else
			{
				mat[i][j] = 0;
			}
		}
	}
	cout << "0000" << endl;
	return ;
}

//******************************************************************************************************************
// Use virtually generated data
//******************************************************************************************************************
void sim_virtual()
{
	char sf[] = "**********************************************************************************************************************************\n";
	char end0[] = "\n***********************************************************************************************************************************\n   ----------------------------------------------------------Sim has done-------------------------------------------------------\n";
	char end1[] = "   ------------------------------------------------------  cycle number is : ";
	char end2[] = "  ---------------------------------------------   \n";
	char end3[] = "   ---------------------------------------------  ";
	char end4[] = "  -------------------------------------------   \n";
	char end5[] = "***********************************************************************************************************************************\n";

	long long cycle = 0;
	long long num_dram_access = 0;
	long long num_sram_access = 0;
	long long prev_cycle = 0;
	long long prev_num_dram_access = 0;
	long long prev_num_sram_access = 0;
	long long compute_latency = 0;
	long long dram_access_latency = 0;

	//initDram();
	//assignDram(DATASET);
	initSram();
	initReg();
	initSim();
	configureWorkload();

	int prev_pcounter_n0 = -1;
	int prev_pcounter_c0 = -1;
	int prev_pcounter_k0 = -1;
	int prev_pcounter_m0 = -1;

	int prev_pcounter_tn0 = -1;
	int prev_pcounter_tc0 = -1;
	int prev_pcounter_tk0 = -1;
	int prev_pcounter_tm0 = -1;

	int prev_pcounter_m1 = -1;
	int prev_pcounter_c1 = -1;
	int prev_pcounter_n1 = -1;
	int prev_pcounter_k1 = -1;

	int prev_pcounter_tm1 = -1;
	int prev_pcounter_tc1 = -1;
	int prev_pcounter_tn1 = -1;
	int prev_pcounter_tk1 = -1;

	clock_t begin_sim, end_sim;
	begin_sim = clock();
	int cycle_Combination = 0;
	int cycle_Aggregation = 0;

	if (!(ACCELERATOR == "GCNAX")) {
		if (ACCELERATOR == "HYGCN") {
			// Notes: the exeuction order of HyGCN is different from that of GCNAX/AWBGCN
			// HyGCN, B[M,K]=A[M,N]*X[N,K], O[M,C]=B[M, K]*W[K, C]; The dimensions of B is different,

			prev_pcounter_m0 = -1;
			prev_pcounter_k0 = -1;
			prev_pcounter_n0 = -1;

			prev_pcounter_tm0 = -1;
			prev_pcounter_tk0 = -1;
			prev_pcounter_tn0 = -1;

			prev_pcounter_m1 = -1;
			prev_pcounter_k1 = -1;
			prev_pcounter_c1 = -1;

			prev_pcounter_tm1 = -1;
			prev_pcounter_tk1 = -1;
			prev_pcounter_tc1 = -1;

			for (pcounter_m0 = 0; pcounter_m0 < DIM_M; pcounter_m0 += TILE_SIZE_TM0) {
				for (pcounter_k0 = 0; pcounter_k0 < DIM_K; pcounter_k0 += TILE_SIZE_TK0) {
					for (pcounter_n0 = 0; pcounter_n0 < DIM_N; pcounter_n0 += TILE_SIZE_TN0) {
						prev_num_dram_access = num_dram_access; //record the current number of dram accesses
						//load SUBMAT_A
						if (pcounter_m0 != prev_pcounter_m0 || pcounter_n0 != prev_pcounter_n0) {
							//load SUBMAT_A
							int num_row = std::min(TILE_SIZE_TM0, DIM_M - pcounter_m0);
							int num_col = std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
							generateVirtualData(SUBMAT_A, num_row, num_col, DENSITY_A);

							//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TK, DIM_K - pcounter_k);
							num_dram_access += countDramAccess(SUBMAT_A);
						}

						//load SUBMAT_W
						if (pcounter_n0 != prev_pcounter_n0 || pcounter_k0 != prev_pcounter_k0) {
							//  load SUBMAT_W
							int num_row = std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
							int num_col = std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0);
							generateVirtualData(SUBMAT_X, num_row, num_col, DENSITY_X);

							//num_dram_access += std::min(TILE_SIZE_TK, DIM_K - pcounter_k) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							num_dram_access += countDramAccess(SUBMAT_X);
						}

						//load SUBMAT_B0
						if (pcounter_m0 != prev_pcounter_m0 || pcounter_k0 != prev_pcounter_k0) {
							//  load SUBMAT_BO
							// if SUBMAT_BO is valid, which means it has to write back to DRAM
							if (valid_submat_bo && prev_pcounter_m0 != -1 && prev_pcounter_k0 != -1) {
								//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
								//num_dram_access += countDramAccess(SUBMAT_X);
								valid_submat_bo = false;
							}

							//  load SUBMAT_BO
							int num_row = std::min(TILE_SIZE_TM0, DIM_M - pcounter_m0);
							int num_col = std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0);
							generateVirtualData(SUBMAT_BO, num_row, num_col, 2);

							//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
						}

						dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access - prev_num_dram_access) /
																   (DRAM_BANDWIDTH /
																	CLOCK_RATE); //calculate the dram access latency

						//on-chip computation (B=AX)
						// loop order tn0->tk->tc0
						for (pcounter_tm0 = 0;
							 pcounter_tm0 < std::min(TILE_SIZE_TM0, DIM_M - pcounter_m0); pcounter_tm0++) {
							for (pcounter_tn0 = 0;
								 pcounter_tn0 < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); pcounter_tn0++) {
								if (SUBMAT_A[pcounter_tm0][pcounter_tn0] != 0) {
									// load element of A to registers
									if (pcounter_tm0 != prev_pcounter_tm0 || pcounter_tn0 != prev_pcounter_tn0) {
										REG_SCALAR = SUBMAT_A[pcounter_tm0][pcounter_tn0];
										num_sram_access += 1;
									}

									for (pcounter_tk0 = 0; pcounter_tk0 < std::min(TILE_SIZE_TK0, DIM_K -
																								  pcounter_k0); pcounter_tk0 += MULTIPLIER_NUM) {
										// load element of X to registers
										if (pcounter_tn0 != prev_pcounter_tn0 || pcounter_tk0 != prev_pcounter_tk0) {
											for (int i = 0;
												 i < std::min(MULTIPLIER_NUM, TILE_SIZE_TK0 - pcounter_tk0); i++) {
												REG_VEC1[i] = SUBMAT_X[pcounter_tn0][pcounter_tk0 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}

										// load elements of B to registers
										if (pcounter_tm0 != prev_pcounter_tm0 || pcounter_tk0 != prev_pcounter_tk0) {
											// if replaced, the results in registers need write back
											if (valid_vec2 && prev_pcounter_tm0 != -1 && prev_pcounter_tk0 != -1) {
												for (int i = 0; i < std::min(MULTIPLIER_NUM,
																			 TILE_SIZE_TK0 - prev_pcounter_tk0); i++) {
													SUBMAT_BO[prev_pcounter_tm0][prev_pcounter_tk0 + i] = REG_VEC2[i];
												}
												valid_vec2 = false;
												num_sram_access += MULTIPLIER_NUM;
											}

											for (int i = 0;
												 i < std::min(MULTIPLIER_NUM, TILE_SIZE_TK0 - pcounter_tk0); i++) {
												REG_VEC2[i] = SUBMAT_BO[pcounter_tm0][pcounter_tk0 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}
										cycle_Aggregation ++;
										compute();
										cycle++;
										//printf(" --------------------------------------------current cycle is---%10d-----------------------------------------------\n", cycle);
										valid_vec2 = true;
										prev_pcounter_tk0 = pcounter_tk0;
										prev_pcounter_tn0 = pcounter_tn0;
										prev_pcounter_tm0 = pcounter_tm0;
									}
								}
							}
						}
						if (valid_vec2 && prev_pcounter_tm0 != -1 && prev_pcounter_tk0 != -1) {
							for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TK0 - prev_pcounter_tk0); i++) {
								SUBMAT_BO[prev_pcounter_tm0][prev_pcounter_tk0 + i] = REG_VEC2[i];
							}
							valid_vec2 = false;
							num_sram_access += MULTIPLIER_NUM;
						}

						compute_latency = cycle - prev_cycle; //recore the latency of computation
						if (dram_access_latency > compute_latency)
							cycle = prev_cycle +
									dram_access_latency; //if dram latency costs more cycles, update current cycle using dram latency

						prev_pcounter_n0 = pcounter_n0;
						valid_submat_bo = true;
						prev_pcounter_k0 = pcounter_k0;
						prev_pcounter_m0 = pcounter_m0;
					}

					for (pcounter_c1 = 0; pcounter_c1 < DIM_C; pcounter_c1 += TILE_SIZE_TC1) {
						prev_num_dram_access = num_dram_access; //record the current number of dram accesses

						// off-chip communication, load SUBMAT_A, SUBMAT_O
						// load SUBMAT_W
						if (pcounter_c1 != prev_pcounter_c1 || pcounter_k0 != prev_pcounter_k0) {
							//  load SUBMAT_A
							int num_row = std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0);
							int num_col = std::min(TILE_SIZE_TC1, DIM_C - pcounter_c1);
							generateVirtualData(SUBMAT_W, num_row, num_col, DENSITY_W);

							//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
							num_dram_access += countDramAccess(SUBMAT_W);
						}

						// load SUBMAT_B0
						if (pcounter_m0 != prev_pcounter_m0 || pcounter_k0 != prev_pcounter_k0) {
							//  load SUBMAT_BI
							int num_row = std::min(TILE_SIZE_TM0, DIM_M - pcounter_m0);
							int num_col = std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0);
							generateVirtualData(SUBMAT_BI, num_row, num_col, 1);

							//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							num_dram_access += countDramAccess(SUBMAT_BI);
						}

						// load SUBMAT_O
						if (pcounter_m0 != prev_pcounter_m0 || pcounter_c1 != prev_pcounter_c1) {
							//  store SUBMAT_O
							// if SUBMAT_BO is valid, which means it has to write back to DRAM
							if (valid_submat_o && prev_pcounter_m0 != -1 && prev_pcounter_c1 != -1) {
								//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
								num_dram_access += countDramAccess(SUBMAT_O);
								valid_submat_o = false;
							}

							//  load SUBMAT_O
							int num_row = std::min(TILE_SIZE_TM0, DIM_M - pcounter_m0);
							int num_col = std::min(TILE_SIZE_TC1, DIM_C - pcounter_c1);
							generateVirtualData(SUBMAT_O, num_row, num_col, 1);

							//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							num_dram_access += countDramAccess(SUBMAT_O);
						}

						dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access - prev_num_dram_access) /
																   (DRAM_BANDWIDTH /
																	CLOCK_RATE); //calculate the dram access latency

						//on-chip computation (O=BW)
						// loop order tn0->tk->tc0
						for (pcounter_tm1 = 0;
							 pcounter_tm1 < std::min(TILE_SIZE_TM0, DIM_M - pcounter_m0); pcounter_tm1++) {
							for (pcounter_tk1 = 0;
								 pcounter_tk1 < std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0); pcounter_tk1++) {
								if (SUBMAT_BO[pcounter_tm1][pcounter_tk1] != 0) {
									// load SUBMAT_B0 to register
									if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tk1 != prev_pcounter_tk1) {
										REG_SCALAR = SUBMAT_BO[pcounter_tm1][pcounter_tk1];
										num_sram_access += 1;
									}

									for (pcounter_tc1 = 0; pcounter_tc1 < std::min(TILE_SIZE_TC1, DIM_C -
																								  pcounter_c1); pcounter_tc1 += MULTIPLIER_NUM) {
										// load element of W to registers
										if (pcounter_tk1 != prev_pcounter_tk1 || pcounter_tc1 != prev_pcounter_tc1) {
											for (int i = 0;
												 i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC1 - pcounter_tc1); i++) {
												REG_VEC1[i] = SUBMAT_W[pcounter_tk1][pcounter_tc1 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}

										// load elements of O to registers
										if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tc1 != prev_pcounter_tc1) {
											// if replaced, the results in registers need write back
											if (valid_vec2 && prev_pcounter_tm1 != -1 && prev_pcounter_tc1 != -1) {
												for (int i = 0; i < std::min(MULTIPLIER_NUM,
																			 TILE_SIZE_TC1 - prev_pcounter_tc1); i++) {
													SUBMAT_O[prev_pcounter_tm1][prev_pcounter_tc1 + i] = REG_VEC2[i];
												}
												valid_vec2 = false;
												num_sram_access += MULTIPLIER_NUM;
											}

											for (int i = 0;
												 i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC1 - pcounter_tc1); i++) {
												REG_VEC2[i] = SUBMAT_O[pcounter_tm1][pcounter_tc1 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}
										cycle_Combination++;
										compute();
										cycle++;
										//printf(" --------------------------------------------current cycle is---%10d-----------------------------------------------\n", cycle);
										valid_vec2 = true;
										prev_pcounter_tc1 = pcounter_tc1;
										prev_pcounter_tk1 = pcounter_tk1;
										prev_pcounter_tm1 = pcounter_tm1;
									}
								}
							}
						}
						if (valid_vec2 && prev_pcounter_tm1 != -1 && prev_pcounter_tc1 != -1) {
							for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC1 - prev_pcounter_tc1); i++) {
								SUBMAT_O[prev_pcounter_tm1][prev_pcounter_tc1 + i] = REG_VEC2[i];
							}
							valid_vec2 = false;
							num_sram_access += MULTIPLIER_NUM;
						}

						compute_latency = cycle - prev_cycle; //recore the latency of computation
						if (dram_access_latency > compute_latency)
							cycle = prev_cycle +
									dram_access_latency; //if dram latency costs more cycles, update current cycle using dram latency

						valid_submat_o = true;
						prev_pcounter_k1 = pcounter_k1;
						prev_pcounter_c1 = pcounter_c1;
						prev_pcounter_m1 = pcounter_m1;
						prev_pcounter_k0 = pcounter_k0;
						prev_pcounter_c0 = pcounter_c0;
						prev_pcounter_m0 = pcounter_m0;
					}
				}
			}
			//  store SUBMAT_O
			// if SUBMAT_O is valid, which means it has to write back to DRAM
			if (valid_submat_o && prev_pcounter_m0 != -1 && prev_pcounter_c1 != -1) {
				//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - prev_pcounter_m) * std::min(TILE_SIZE_TC0, DIM_C - prev_pcounter_c1);
				num_dram_access += countDramAccess(SUBMAT_O);
			}
			valid_submat_o = false;
			cout<<"Agg:"<<cycle_Aggregation<<endl;
			cout<<"Comb:"<<cycle_Combination<<endl;
		} else if (ACCELERATOR == "AWBGCN") {
			COMPUTE_TYPE = 0;

			prev_pcounter_n0 = -1;
			prev_pcounter_c0 = -1;
			prev_pcounter_k0 = -1;
			prev_pcounter_tn0 = -1;
			prev_pcounter_tc0 = -1;
			prev_pcounter_tk0 = -1;

			prev_pcounter_m1 = -1;
			prev_pcounter_c1 = -1;
			prev_pcounter_n1 = -1;
			prev_pcounter_tm1 = -1;
			prev_pcounter_tc1 = -1;
			prev_pcounter_tn1 = -1;

			for (pcounter_n0 = 0; pcounter_n0 < DIM_N; pcounter_n0 += TILE_SIZE_TN0) {
				for (pcounter_c0 = 0; pcounter_c0 < DIM_C; pcounter_c0 += TILE_SIZE_TC0) {
					for (pcounter_k0 = 0; pcounter_k0 < DIM_K; pcounter_k0 += TILE_SIZE_TK0) {
						prev_num_dram_access = num_dram_access; //record the current number of dram accesses
						if (pcounter_n0 != prev_pcounter_n0 || pcounter_k0 != prev_pcounter_k0) {
							//  load SUBMAT_X
							int num_row = std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
							int num_col = std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0);
							generateVirtualData(SUBMAT_X, num_row, num_col, DENSITY_X);

							//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TK, DIM_K - pcounter_k);
							num_dram_access += countDramAccess(SUBMAT_X);
						}

						if (pcounter_k0 != prev_pcounter_k0 || pcounter_c0 != prev_pcounter_c0) {
							//  load SUBMAT_W
							int num_row = std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0);
							int num_col = std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							generateVirtualData(SUBMAT_W, num_row, num_col, DENSITY_W);

							//num_dram_access += std::min(TILE_SIZE_TK, DIM_K - pcounter_k) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							num_dram_access += countDramAccess(SUBMAT_W);
						}

						if (pcounter_n0 != prev_pcounter_n0 || pcounter_c0 != prev_pcounter_c0) {
							//  load SUBMAT_BO
							// if SUBMAT_BO is valid, which means it has to write back to DRAM
							if (valid_submat_bo && prev_pcounter_n0 != -1 && prev_pcounter_c0 != -1) {
								//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
								//num_dram_access += countDramAccess(SUBMAT_X);
								valid_submat_bo = false;
							}

							//  load SUBMAT_BO
							int num_row = std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
							int num_col = std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							generateVirtualData(SUBMAT_BO, num_row, num_col, 1);

							//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
						}

						dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access - prev_num_dram_access) /
																   (DRAM_BANDWIDTH /
																	CLOCK_RATE); //calculate the dram access latency

						//on-chip computation
						// loop order tn0->tk->tc0
						for (pcounter_tn0 = 0;
							 pcounter_tn0 < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); pcounter_tn0++) {
							for (pcounter_tc0 = 0;
								 pcounter_tc0 < std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0); pcounter_tc0++) {
								if (true) // use inner product, cannot skip zero computations
								{

									if (pcounter_tn0 != prev_pcounter_tn0 || pcounter_tc0 != prev_pcounter_tc0) {
										// if the element has to be written back to SRAM
										if (valid_scalar && prev_pcounter_tn0 != -1 && prev_pcounter_tc0 != -1) {
											SUBMAT_BO[prev_pcounter_tn0][prev_pcounter_tc0] = REG_SCALAR;
										}

										REG_SCALAR = SUBMAT_BO[pcounter_tn0][pcounter_tc0];
										num_sram_access += 1;
									}

									for (pcounter_tk0 = 0; pcounter_tk0 < std::min(TILE_SIZE_TK0, DIM_K -
																								  pcounter_k0); pcounter_tk0 += MULTIPLIER_NUM) {
										// load element of W to registers
										if (pcounter_tk0 != prev_pcounter_tk0 || pcounter_tc0 != prev_pcounter_tc0) {
											for (int i = 0;
												 i < std::min(MULTIPLIER_NUM, TILE_SIZE_TK0 - pcounter_tk0); i++) {
												REG_VEC1[i] = SUBMAT_W[pcounter_tk0 + i][pcounter_tc0];
											}
											num_sram_access += MULTIPLIER_NUM;
										}

										// load elements of B to registers
										if (pcounter_tn0 != prev_pcounter_tn0 || pcounter_tk0 != prev_pcounter_tk0) {

											for (int i = 0;
												 i < std::min(MULTIPLIER_NUM, TILE_SIZE_TK0 - pcounter_tk0); i++) {
												REG_VEC2[i] = SUBMAT_X[pcounter_tn0][pcounter_tk0 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}

										int num_nonzero_computation = 0;
										for (int i = 0; i < MULTIPLIER_NUM; i++) {
											if (REG_VEC1[i] != 0 || REG_VEC2[i] != 0)
												num_nonzero_computation++;
										}

										compute();
										cycle += num_nonzero_computation;
										//printf(" --------------------------------------------current cycle is---%10d-----------------------------------------------\n", cycle);
										valid_vec2 = true;
										prev_pcounter_tc0 = pcounter_tc0;
										prev_pcounter_tk0 = pcounter_tk0;
										prev_pcounter_tn0 = pcounter_tn0;
									}
								}
							}
						}
						if (valid_scalar && prev_pcounter_tn0 != -1 && prev_pcounter_tc0 != -1) {
							SUBMAT_BO[prev_pcounter_tn0][prev_pcounter_tc0] = REG_SCALAR;
							valid_scalar = false;
							num_sram_access += 1;
						}

						compute_latency = cycle - prev_cycle; //recore the latency of computation
						if (dram_access_latency > compute_latency / MULTIPLIER_NUM)
							cycle = prev_cycle + dram_access_latency *
												 MULTIPLIER_NUM; //if dram latency costs more cycles, update current cycle using dram latency

						prev_pcounter_k0 = pcounter_k0;
						valid_submat_bo = true;
						prev_pcounter_c0 = pcounter_c0;
						prev_pcounter_n0 = pcounter_n0;
					}

					for (pcounter_m1 = 0; pcounter_m1 < DIM_M; pcounter_m1 += TILE_SIZE_TM1) {
						prev_num_dram_access = num_dram_access; //record the current number of dram accesses

						// off-chip communication, load SUBMAT_A, SUBMAT_O
						if (pcounter_m1 != prev_pcounter_m1 || pcounter_n0 != prev_pcounter_n0) {
							//  load SUBMAT_A
							int num_row = std::min(TILE_SIZE_TM0, DIM_M - pcounter_m1);
							int num_col = std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
							generateVirtualData(SUBMAT_A, num_row, num_col, DENSITY_A);

							//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
							num_dram_access += countDramAccess(SUBMAT_A);
						}

						if (pcounter_n0 != prev_pcounter_n0 || pcounter_c0 != prev_pcounter_c0) {
							//  load SUBMAT_BI
							int num_row = std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
							int num_col = std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							generateVirtualData(SUBMAT_BI, num_row, num_col, 1);

							//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							num_dram_access += countDramAccess(SUBMAT_BI);
						}

						if (pcounter_m1 != prev_pcounter_m1 || pcounter_c0 != prev_pcounter_c0) {
							//  store SUBMAT_O
							// if SUBMAT_BO is valid, which means it has to write back to DRAM
							if (valid_submat_o && prev_pcounter_m1 != -1 && prev_pcounter_c0 != -1) {
								//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
								num_dram_access += countDramAccess(SUBMAT_O);
								valid_submat_o = false;
							}

							//  load SUBMAT_O
							int num_row = std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1);
							int num_col = std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							generateVirtualData(SUBMAT_O, num_row, num_col, 1);

							//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							num_dram_access += countDramAccess(SUBMAT_O);
						}

						dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access - prev_num_dram_access) /
																   (DRAM_BANDWIDTH /
																	CLOCK_RATE); //calculate the dram access latency

						//on-chip computation
						// loop order tn0->tk->tc0
						for (pcounter_tm1 = 0;
							 pcounter_tm1 < std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1); pcounter_tm1++) {
							for (pcounter_tc1 = 0;
								 pcounter_tc1 < std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0); pcounter_tc1++) {
								if (true) {
									if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tc1 != prev_pcounter_tc1) {
										// if the element has to be written back to SRAM
										if (valid_scalar && prev_pcounter_tm1 != -1 && prev_pcounter_tc1 != -1) {
											SUBMAT_O[prev_pcounter_tm1][prev_pcounter_tc1] = REG_SCALAR;
										}
										REG_SCALAR = SUBMAT_O[pcounter_tm1][pcounter_tc1];
										num_sram_access += 1;
									}

									for (pcounter_tn1 = 0; pcounter_tn1 < std::min(TILE_SIZE_TN0, DIM_C -
																								  pcounter_n0); pcounter_tn1 += MULTIPLIER_NUM) {
										// load element of BI to registers
										if (pcounter_tn1 != prev_pcounter_tn1 || pcounter_tc1 != prev_pcounter_tc1) {
											for (int i = 0;
												 i < std::min(MULTIPLIER_NUM, TILE_SIZE_TN0 - pcounter_tn1); i++) {
												REG_VEC1[i] = SUBMAT_BO[pcounter_tn1 + i][pcounter_tc1];
											}
											num_sram_access += MULTIPLIER_NUM;
										}

										// load elements of O to registers
										if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tn1 != prev_pcounter_tn1) {

											for (int i = 0;
												 i < std::min(MULTIPLIER_NUM, TILE_SIZE_TN0 - pcounter_tn1); i++) {
												REG_VEC2[i] = SUBMAT_A[pcounter_tm1][pcounter_tn1 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}
										int num_nonzero_computation = 0;
										for (int i = 0; i < MULTIPLIER_NUM; i++) {
											if (REG_VEC1[i] != 0 || REG_VEC2[i] != 0)
												num_nonzero_computation++;
										}

										compute();
										cycle += num_nonzero_computation;

										//printf(" --------------------------------------------current cycle is---%10d-----------------------------------------------\n", cycle);
										valid_vec2 = true;
										prev_pcounter_tc1 = pcounter_tc1;
										prev_pcounter_tn1 = pcounter_tn1;
										prev_pcounter_tm1 = pcounter_tm1;
									}
								}
							}
						}
						if (valid_scalar && prev_pcounter_tm1 != -1 && prev_pcounter_tc1 != -1) {
							SUBMAT_O[prev_pcounter_tm1][prev_pcounter_tc1] = REG_SCALAR;
							valid_scalar = false;
							num_sram_access += 1;
						}

						compute_latency = cycle - prev_cycle; //recore the latency of computation
						if (dram_access_latency > compute_latency / MULTIPLIER_NUM)
							cycle = prev_cycle + dram_access_latency *
												 MULTIPLIER_NUM; //if dram latency costs more cycles, update current cycle using dram latency

						valid_submat_o = true;
						prev_pcounter_n1 = pcounter_n1;
						prev_pcounter_c1 = pcounter_c1;
						prev_pcounter_m1 = pcounter_m1;
						prev_pcounter_n0 = pcounter_n0;
						prev_pcounter_c0 = pcounter_c0;
						prev_pcounter_m0 = pcounter_m0;
					}
				}
			}
			//  store SUBMAT_O
			// if SUBMAT_O is valid, which means it has to write back to DRAM
			if (valid_submat_o && prev_pcounter_m1 != -1 && prev_pcounter_c0 != -1) {
				//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - prev_pcounter_m) * std::min(TILE_SIZE_TC0, DIM_C - prev_pcounter_c1);
				num_dram_access += countDramAccess(SUBMAT_O);
			}
			valid_submat_o = false;
			cycle = cycle / MULTIPLIER_NUM;
		} else if (ACCELERATOR == "HDGCN") {
			// personal used parameter
			int fold = 0;
			int fold_DramAccessofWeight = 0;
			int flag = 0;
			int tmp = 0;
			int SA_row = sqrt(MULTIPLIER_NUM);
			int SA_col = sqrt(MULTIPLIER_NUM);

			prev_pcounter_n0 = -1;
			prev_pcounter_c0 = -1;
			prev_pcounter_k0 = -1;
			prev_pcounter_tn0 = -1;
			prev_pcounter_tc0 = -1;
			prev_pcounter_tk0 = -1;

			prev_pcounter_m1 = -1;
			prev_pcounter_c1 = -1;
			prev_pcounter_n1 = -1;
			prev_pcounter_tm1 = -1;
			prev_pcounter_tc1 = -1;
			prev_pcounter_tn1 = -1;

			for (pcounter_n0 = 0; pcounter_n0 < DIM_N; pcounter_n0 += TILE_SIZE_TN0) {
				for (pcounter_c0 = 0; pcounter_c0 < DIM_C; pcounter_c0 += TILE_SIZE_TC0) {
					// Combination Phase: B=XW
					for (pcounter_k0 = 0; pcounter_k0 < DIM_K; pcounter_k0 += TILE_SIZE_TK0) {
						prev_num_dram_access = num_dram_access; //record the current number of dram accesses

						//  load SUBMAT_X
						if (pcounter_n0 != prev_pcounter_n0 || pcounter_k0 != prev_pcounter_k0) {
							
							int num_row = std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
							int num_col = std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0);
							generateVirtualData(SUBMAT_X, num_row, num_col, DENSITY_X);

							num_dram_access += countDramAccess(SUBMAT_X);
						}

						//  load SUBMAT_W
						if (pcounter_k0 != prev_pcounter_k0 || pcounter_c0 != prev_pcounter_c0) {
							
							int num_row = std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0);
							int num_col = std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							generateVirtualData(SUBMAT_W, num_row, num_col, DENSITY_W);

							//num_dram_access += std::min(TILE_SIZE_TK, DIM_K - pcounter_k) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							tmp = countDramAccess(SUBMAT_W);
							num_dram_access += tmp;
							fold_DramAccessofWeight += tmp;
						}

						//  load SUBMAT_B
						if (pcounter_n0 != prev_pcounter_n0 || pcounter_c0 != prev_pcounter_c0) {
							//  load SUBMAT_BO
							// if SUBMAT_BO is valid, which means it has to write back to DRAM
							if (valid_submat_bo && prev_pcounter_n0 != -1 && prev_pcounter_c0 != -1) {
								//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
								//num_dram_access += countDramAccess(SUBMAT_X);
								valid_submat_bo = false;
							}

							//  load SUBMAT_BO
							int num_row = std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
							int num_col = std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							generateVirtualData(SUBMAT_BO, num_row, num_col, 1);

							//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
						}

						dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access - prev_num_dram_access) /
																   (DRAM_BANDWIDTH /
																	CLOCK_RATE); //calculate the dram access latency

						// used to calculate the redundant sram_access of weight
						int fold_weight_sram_access = 0;
						int fold_weight = 0;
						int sram_access_weight = 0;

						//on-chip computation
						// loop order tn0->tk->tc0
						for (pcounter_tn0 = 0;
							 pcounter_tn0 < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); pcounter_tn0++) {
							for (pcounter_tk0 = 0;
								 pcounter_tk0 < std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0); pcounter_tk0++) {
								// load element of X to registers
								if (SUBMAT_X[pcounter_tn0][pcounter_tk0] != 0) {
									if (pcounter_tn0 != prev_pcounter_tn0 || pcounter_tk0 != prev_pcounter_tk0) {
										REG_SCALAR = SUBMAT_X[pcounter_tn0][pcounter_tk0];
										num_sram_access += 1;
										//num_sram_access += sqrt(MULTIPLIER_NUM);
									}

									for (pcounter_tc0 = 0; pcounter_tc0 < std::min(TILE_SIZE_TC0, DIM_C -
																								  pcounter_c0); pcounter_tc0 += MULTIPLIER_NUM) {
										// load element of W to registers
										if (pcounter_tk0 != prev_pcounter_tk0 || pcounter_tc0 != prev_pcounter_tc0) {
											for (int i = 0;
												 i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc0); i++) {
												REG_VEC1[i] = SUBMAT_W[pcounter_tk0][pcounter_tc0 + i];
											}
											// used to record the reuse of weight matrix
											num_sram_access += MULTIPLIER_NUM;
											sram_access_weight += MULTIPLIER_NUM;
										}

										// load elements of B to registers
										if (pcounter_tn0 != prev_pcounter_tn0 || pcounter_tc0 != prev_pcounter_tc0) {
											// if replaced, the results in registers need write back
											if (valid_vec2 && prev_pcounter_tn0 != -1 && prev_pcounter_tc0 != -1) {
												for (int i = 0; i < std::min(MULTIPLIER_NUM,
																			 TILE_SIZE_TC0 - prev_pcounter_tc0); i++) {
													SUBMAT_BO[prev_pcounter_tn0][prev_pcounter_tc0 + i] = REG_VEC2[i];
												}
												valid_vec2 = false;
												num_sram_access += MULTIPLIER_NUM;
											}

											for (int i = 0;
												 i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc0); i++) {
												REG_VEC2[i] = SUBMAT_BO[pcounter_tn0][pcounter_tc0 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}
										compute();
										cycle++;
										cycle_Combination++;
										//printf(" --------------------------------------------current cycle is---%10d-----------------------------------------------\n", cycle);
										valid_vec2 = true;
										prev_pcounter_tc0 = pcounter_tc0;
										prev_pcounter_tk0 = pcounter_tk0;
										prev_pcounter_tn0 = pcounter_tn0;
									}
								}
							}
							fold_weight++;
							fold_weight_sram_access = sram_access_weight;
							sram_access_weight = 0;
						}

						// sram reduction because of weight reuse
						num_sram_access -= (fold_weight-1)*fold_weight_sram_access;
					
						if (valid_vec2 && prev_pcounter_tn0 != -1 && prev_pcounter_tc0 != -1) {
							for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - prev_pcounter_tc0); i++) {
								SUBMAT_BO[prev_pcounter_tn0][prev_pcounter_tc0 + i] = REG_VEC2[i];
							}
							valid_vec2 = false;
							num_sram_access += MULTIPLIER_NUM;
						}

						compute_latency = cycle - prev_cycle; //recore the latency of computation
						if (dram_access_latency > compute_latency)
							cycle = prev_cycle +
									dram_access_latency; //if dram latency costs more cycles, update current cycle using dram latency

						prev_pcounter_k0 = pcounter_k0;
						valid_submat_bo = true;
						prev_pcounter_c0 = pcounter_c0;
						prev_pcounter_n0 = pcounter_n0;
					}

					
					// Aggregation Phase: O=AB
					for (pcounter_m1 = 0; pcounter_m1 < DIM_M; pcounter_m1 += TILE_SIZE_TM1) {
						prev_num_dram_access = num_dram_access; //record the current number of dram accesses

						// off-chip communication, load SUBMAT_A, SUBMAT_O
						// load SUBMAT_A
						if (pcounter_m1 != prev_pcounter_m1 || pcounter_n0 != prev_pcounter_n0) {
							//  load SUBMAT_A
							int num_row = std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1);
							int num_col = std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
							generateVirtualData(SUBMAT_A, num_row, num_col, DENSITY_A);

							//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
							num_dram_access += countDramAccess(SUBMAT_A);
						}

						// load SUBMAT_B
						if (pcounter_n0 != prev_pcounter_n0 || pcounter_c0 != prev_pcounter_c0) {
							//  load SUBMAT_BI
							int num_row = std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
							int num_col = std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							generateVirtualData(SUBMAT_BI, num_row, num_col, 1);

							//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							//num_dram_access += countDramAccess(SUBMAT_BI);
						}

						// load SUBMAT_O
						if (pcounter_m1 != prev_pcounter_m1 || pcounter_c0 != prev_pcounter_c0) {
							//  store SUBMAT_O
							// if SUBMAT_BO is valid, which means it has to write back to DRAM
							if (valid_submat_o && prev_pcounter_m1 != -1 && prev_pcounter_c0 != -1) {
								//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
								num_dram_access += countDramAccess(SUBMAT_O);
								valid_submat_o = false;
							}

							//  load SUBMAT_O
							int num_row = std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1);
							int num_col = std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							generateVirtualData(SUBMAT_O, num_row, num_col, 1);

							//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							num_dram_access += countDramAccess(SUBMAT_O);
						}

						dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access - prev_num_dram_access) /
																   (DRAM_BANDWIDTH /
																	CLOCK_RATE); //calculate the dram access latency

						int sram_B_access = 0;
						int fold_sram_B_access = 0;
						int fold_B_sram_access = 0;


						//on-chip computation
						// loop order tn0->tk->tc0
						for (pcounter_tm1 = 0;
							 pcounter_tm1 < std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1); pcounter_tm1++) {
							for (pcounter_tn1 = 0;
								 pcounter_tn1 < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); pcounter_tn1++) {
								if (SUBMAT_A[pcounter_tm1][pcounter_tn1] != 0) {
									// load element of A to registers
									if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tn1 != prev_pcounter_tn1) {
										REG_SCALAR = SUBMAT_A[pcounter_tm1][pcounter_tn1];
										num_sram_access += 1;
									}

									for (pcounter_tc1 = 0; pcounter_tc1 < std::min(TILE_SIZE_TC0, DIM_C -
																								  pcounter_c0); pcounter_tc1 += MULTIPLIER_NUM) {
										// load element of BI to registers
										if (pcounter_tn1 != prev_pcounter_tn1 || pcounter_tc1 != prev_pcounter_tc1) {
											for (int i = 0;
												 i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc1); i++) {
												REG_VEC1[i] = SUBMAT_BO[pcounter_tn1][pcounter_tc1 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
											sram_B_access += MULTIPLIER_NUM;
										}

										// load elements of O to registers
										if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tc1 != prev_pcounter_tc1) {
											// if replaced, the results in registers need write back
											if (valid_vec2 && prev_pcounter_tm1 != -1 && prev_pcounter_tc1 != -1) {
												for (int i = 0; i < std::min(MULTIPLIER_NUM,
																			 TILE_SIZE_TC0 - prev_pcounter_tc1); i++) {
													SUBMAT_O[prev_pcounter_tm1][prev_pcounter_tc1 + i] = REG_VEC2[i];
												}
												valid_vec2 = false;
												num_sram_access += MULTIPLIER_NUM;
											}

											for (int i = 0;
												 i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc1); i++) {
												REG_VEC2[i] = SUBMAT_O[pcounter_tm1][pcounter_tc1 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}
										compute();
										cycle_Aggregation++;
										cycle++;
										//printf(" --------------------------------------------current cycle is---%10d-----------------------------------------------\n", cycle);
										valid_vec2 = true;
										prev_pcounter_tc1 = pcounter_tc1;
										prev_pcounter_tn1 = pcounter_tn1;
										prev_pcounter_tm1 = pcounter_tm1;
									}
								}
							}
							fold_B_sram_access++;
							fold_B_sram_access = sram_B_access;
							sram_B_access = 0;
						}

						// sram access reduction because of reuse of Matrix B
						// num_sram_access -= fold_B_sram_access * (fold_B_sram_access-1);

						if (valid_vec2 && prev_pcounter_tm1 != -1 && prev_pcounter_tc1 != -1) {
							for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - prev_pcounter_tc1); i++) {
								SUBMAT_O[prev_pcounter_tm1][prev_pcounter_tc1 + i] = REG_VEC2[i];
							}
							valid_vec2 = false;
							num_sram_access += MULTIPLIER_NUM;
						}

						compute_latency = cycle - prev_cycle; //recore the latency of computation
						if (dram_access_latency > compute_latency)
							cycle = prev_cycle +
									dram_access_latency; //if dram latency costs more cycles, update current cycle using dram latency

						valid_submat_o = true;
						prev_pcounter_n1 = pcounter_n1;
						prev_pcounter_c1 = pcounter_c1;
						prev_pcounter_m1 = pcounter_m1;
						prev_pcounter_n0 = pcounter_n0;
						prev_pcounter_c0 = pcounter_c0;
						prev_pcounter_m0 = pcounter_m0;
					}
				}

				if (pcounter_n0 == 0) {  flag = fold_DramAccessofWeight; }
				fold_DramAccessofWeight = 0;
				fold++;
			}

			// Dram access reduction because of weight reuse
			cout << fold << endl;
			num_dram_access = num_dram_access - (fold-1) * flag;

			//  store SUBMAT_O
			// if SUBMAT_O is valid, which means it has to write back to DRAM
			if (valid_submat_o && prev_pcounter_m1 != -1 && prev_pcounter_c0 != -1) {
				//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - prev_pcounter_m) * std::min(TILE_SIZE_TC0, DIM_C - prev_pcounter_c1);
				num_dram_access += countDramAccess(SUBMAT_O);
			}
			valid_submat_o = false;
			cout<<"comb:"<<ceil(cycle_Combination * 11 / 13)<<endl;
			cout<<"AGG:"<<ceil(cycle_Aggregation * 11 / 13) <<endl;

			cycle = ceil( cycle * 11 / 13 );

		} else if (ACCELERATOR == "EGCNAX") {
			// personal parameters
			prev_pcounter_n0 = -1;
            prev_pcounter_c0 = -1;
            prev_pcounter_k0 = -1;
            prev_pcounter_tn0 = -1;
            prev_pcounter_tc0 = -1;
            prev_pcounter_tk0 = -1;

            prev_pcounter_m1 = -1;
            prev_pcounter_c1 = -1;
            prev_pcounter_n1 = -1;
            prev_pcounter_tm1 = -1;
            prev_pcounter_tc1 = -1;
            prev_pcounter_tn1 = -1;

			// check perform the aggregation first or combination
			if (DIM_K > DIM_C) {
				// perform the combination first and then aggregation

				// traverse the CSR vector and sort the vertices from hot vertices to cold vertices
				cycle += DIM_N;

				// record the DRAM access of loading the related data from DRAM to the on-chip buffer
				int pre_row = int (0.2 * DIM_N);
				int pre_col = DIM_K;
				cout << pre_row << " " << pre_col << endl;
				cout << "11" << endl;
				generateVirtualData(SUBMAT_X, pre_row, pre_col, DENSITY_X);
				cout << "222" << endl;
				num_dram_access += countDramAccess(SUBMAT_X);
				return;

				// record the new density of the A matrix based on the power-law
				double den_temp_A = 0.0;
				int a = DIM_M * DIM_N * DENSITY_A;
				a = 0.2 * a;
				den_temp_A = a / (0.8 * DIM_M * DIM_N);
				
				// time overhead of the pre-loading data
				dram_access_latency = DRAM_SETUP_LATENCY + (pre_col * pre_row) / (DRAM_BANDWIDTH / CLOCK_RATE);
				cycle += dram_access_latency;

				// off-chip communication: B=X*W
				// load matrix X, W, B to SUBMAT_X, SUBMAT_W, SUBMAT_BO
				for (pcounter_n0 = 0; pcounter_n0 < DIM_N; pcounter_n0 += TILE_SIZE_TN0) {
					for (pcounter_c0 = 0; pcounter_c0 < DIM_C; pcounter_c0 += TILE_SIZE_TC0) {
						for (pcounter_k0 = 0; pcounter_k0 < DIM_K; pcounter_k0 += TILE_SIZE_TK0) {

							// record the current number of dram access
							prev_num_dram_access = num_dram_access;

							// load SUBMAT_X
							if (pcounter_n0 < 0.2 * DIM_N) {
								// hot vertices are already pre-loaded
								if (pcounter_n0 != prev_pcounter_n0 || pcounter_k0 != prev_pcounter_k0) {
									int num_row = std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
									int num_col = std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0);
									generateVirtualData(SUBMAT_X, num_row, num_col, DENSITY_X);
								}
							} else {
								// vertices are not pre-loaded
								if (pcounter_n0 != prev_pcounter_n0 || pcounter_k0 != prev_pcounter_k0) {
									int num_row = std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
									int num_col = std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0);
									generateVirtualData(SUBMAT_X, num_row, num_col, DENSITY_X);

									num_dram_access += countDramAccess(SUBMAT_X);
								}
							}

							// load SUBMAT_W
							if (pcounter_k0 != prev_pcounter_k0 || pcounter_c0 != prev_pcounter_c0) {
								int num_row = std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0);
								int num_col = std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
								generateVirtualData(SUBMAT_W, num_row, num_col, DENSITY_W);

								num_dram_access += countDramAccess(SUBMAT_W);
							}

							// load SUBMAT_BO
							if (pcounter_n0 != prev_pcounter_n0 || pcounter_c0 != prev_pcounter_c0) {
								// if SUBMAT_BO is valid, which means it has to write back to DRAM
								if (valid_submat_bo && prev_pcounter_n0 != -1 && prev_pcounter_c0 != -1) {
									if (pcounter_n0 < 0.2 * DIM_N) {
										// write back to on-chip buffer
										num_sram_access += countDramAccess(SUBMAT_BO);
										valid_submat_bo = false;
									} else {
										// stream SUBMAT_BO back to DRAM
										num_dram_access += countDramAccess(SUBMAT_BO);
										valid_submat_bo = false;
									}
									
								}

								// if SUBMAT_BO is invalid, load data from the DRAM
								if (pcounter_n0 < 0.2 * DIM_N) {
									// the needed data are already pre-loaded
									int num_row = std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0);
									int num_col = std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
									generateVirtualData(SUBMAT_BO, num_row, num_col, 1);
								} else {
									int num_row = std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0);
									int num_col = std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
									generateVirtualData(SUBMAT_BO, num_row, num_col, 1);

									num_dram_access += countDramAccess(SUBMAT_BO);
								}
							}

							// calculate the dram access latency
							dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access-prev_num_dram_access) / (DRAM_BANDWIDTH/CLOCK_RATE);

							//systolic array uses WS for processing matrix multiplication
							//number of weight sram access
							for (pcounter_tn0 = 0; pcounter_tn0 < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); pcounter_tn0++) {
								for (pcounter_tk0 = 0; pcounter_tk0 < std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0); pcounter_tk0++) {
									if (SUBMAT_X[pcounter_tn0][pcounter_tk0] != 0) {
										if (pcounter_tn0 != prev_pcounter_tn0 || pcounter_tk0 != prev_pcounter_tk0) {
											REG_SCALAR = SUBMAT_X[pcounter_tn0][pcounter_tk0];
											num_sram_access += 1;
										}

										for (pcounter_tc0 = 0; pcounter_tc0 < std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0); pcounter_tc0 += MULTIPLIER_NUM) {
											// load element of W to registers
											if (pcounter_tk0 != prev_pcounter_tk0 || pcounter_tc0 != prev_pcounter_tc0) {
												for (int i = 0;
													i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc0); i++) {
													REG_VEC1[i] = SUBMAT_W[pcounter_tk0][pcounter_tc0 + i];
												}
												num_sram_access += MULTIPLIER_NUM;
											}

											// load elements of B to registers
											if (pcounter_tn0 != prev_pcounter_tn0 || pcounter_tc0 != prev_pcounter_tc0) {
												// if replaced, the results in registers need write back
												if (valid_vec2 && prev_pcounter_tn0 != -1 && prev_pcounter_tc0 != -1) {
													for (int i = 0; i < std::min(MULTIPLIER_NUM,
																				TILE_SIZE_TC0 - prev_pcounter_tc0); i++) {
														SUBMAT_BO[prev_pcounter_tn0][prev_pcounter_tc0 + i] = REG_VEC2[i];
													}
													valid_vec2 = false;
													num_sram_access += MULTIPLIER_NUM;
												}

												for (int i = 0;
													i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc0); i++) {
													REG_VEC2[i] = SUBMAT_BO[pcounter_tn0][pcounter_tc0 + i];
												}
												num_sram_access += MULTIPLIER_NUM;
											}
											compute();
											cycle = cycle + 0.9;
											//printf(" --------------------------------------------current cycle is---%10d-----------------------------------------------\n", cycle);
											valid_vec2 = true;
											prev_pcounter_tc0 = pcounter_tc0;
											prev_pcounter_tk0 = pcounter_tk0;
											prev_pcounter_tn0 = pcounter_tn0;
										}
									}
								}
							}

							if (valid_vec2 && prev_pcounter_tn0 != -1 && prev_pcounter_tc0 != -1) {
								for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - prev_pcounter_tc0); i++) {
									SUBMAT_BO[prev_pcounter_tn0][prev_pcounter_tc0 + i] = REG_VEC2[i];
								}
								valid_vec2 = false;
								num_sram_access += MULTIPLIER_NUM;
							}

							compute_latency = cycle - prev_cycle;
							if (dram_access_latency > compute_latency) {
								cycle = prev_cycle + dram_access_latency;
							}
							valid_submat_bo = true;
							prev_pcounter_k0 = pcounter_k0;
							prev_pcounter_c0 = pcounter_c0;
							prev_pcounter_n0 = pcounter_n0;
						}

						cout << num_dram_access << endl;
						for (pcounter_m1 = 0; pcounter_m1 < DIM_M; pcounter_m1 += TILE_SIZE_TM1) {
							prev_num_dram_access = num_dram_access;

							// load SUBMAT_A
							if (pcounter_m1 != prev_pcounter_m1 || pcounter_n0 != prev_pcounter_n0) {
								int num_row = std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1);
								int num_col = std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
								generateVirtualData(SUBMAT_A, num_row, num_col, DENSITY_A);

								//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
								num_dram_access += countDramAccess(SUBMAT_A);
							}

							// load SUBMAT_BI
							if (pcounter_n0 != prev_pcounter_n0 || pcounter_c0 != prev_pcounter_c0) {
								// check if the data is pre-loaded
								if (pcounter_n0 < 0.2 * DIM_N) {
									int num_row = std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
									int num_col = std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
									generateVirtualData(SUBMAT_BI, num_row, num_col, 1);
								} else {
									int num_row = std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
									int num_col = std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
									generateVirtualData(SUBMAT_BI, num_row, num_col, 1);

									num_dram_access += countDramAccess(SUBMAT_BI);
								}
							}

							// load SUBMAT_O
							if (pcounter_m1 != prev_pcounter_m1 || pcounter_c0 != prev_pcounter_c0) {
								// check if the data needs to write back
								if (valid_submat_o && prev_pcounter_m1 != -1 && prev_pcounter_c0 != -1) {
									//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
									num_dram_access += countDramAccess(SUBMAT_O);
									valid_submat_o = false;
								}

								// check if the data is pre-loaded
								if (pcounter_m1 < 0.2 * DIM_M) {
									int num_row = std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1);
									int num_col = std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
									generateVirtualData(SUBMAT_O, num_row, num_col, 1);
								} else {
									int num_row = std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1);
									int num_col = std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
									generateVirtualData(SUBMAT_O, num_row, num_col, 1);

									num_dram_access += countDramAccess(SUBMAT_O);
								}
							}

							dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access - prev_num_dram_access) / (DRAM_BANDWIDTH / CLOCK_RATE); //calculate the dram access latency

							// on-chip computation
							for (pcounter_tm1 = 0;
								pcounter_tm1 < std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1); pcounter_tm1++) {
								for (pcounter_tn1 = 0;
									pcounter_tn1 < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); pcounter_tn1++) {
									if (SUBMAT_A[pcounter_tm1][pcounter_tn1] != 0) {
										if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tn1 != prev_pcounter_tn1) {
											REG_SCALAR = SUBMAT_A[pcounter_tm1][pcounter_tn1];
											num_sram_access += 1;
										}

										for (pcounter_tc1 = 0; pcounter_tc1 < std::min(TILE_SIZE_TC0, DIM_C -
																									pcounter_c0); pcounter_tc1 += MULTIPLIER_NUM) {
											// load element of BI to registers
											if (pcounter_tn1 != prev_pcounter_tn1 || pcounter_tc1 != prev_pcounter_tc1) {
												for (int i = 0;
													i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc1); i++) {
													REG_VEC1[i] = SUBMAT_BO[pcounter_tn1][pcounter_tc1 + i];
												}
												num_sram_access += MULTIPLIER_NUM;
											}

											// load elements of O to registers
											if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tc1 != prev_pcounter_tc1) {
												// if replaced, the results in registers need write back
												if (valid_vec2 && prev_pcounter_tm1 != -1 && prev_pcounter_tc1 != -1) {
													for (int i = 0; i < std::min(MULTIPLIER_NUM,
																				TILE_SIZE_TC0 - prev_pcounter_tc1); i++) {
														SUBMAT_O[prev_pcounter_tm1][prev_pcounter_tc1 + i] = REG_VEC2[i];
													}
													valid_vec2 = false;
													num_sram_access += MULTIPLIER_NUM;
												}

												for (int i = 0;
													i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc1); i++) {
													REG_VEC2[i] = SUBMAT_O[pcounter_tm1][pcounter_tc1 + i];
												}
												num_sram_access += MULTIPLIER_NUM;
											}
											compute();
											cycle++;
											//printf(" --------------------------------------------current cycle is---%10d-----------------------------------------------\n", cycle);
											valid_vec2 = true;
											prev_pcounter_tc1 = pcounter_tc1;
											prev_pcounter_tn1 = pcounter_tn1;
											prev_pcounter_tm1 = pcounter_tm1;
										}
									}
								}
							}

							if (valid_vec2 && prev_pcounter_tm1 != -1 && prev_pcounter_tc1 != -1) {
								for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - prev_pcounter_tc1); i++) {
									SUBMAT_O[prev_pcounter_tm1][prev_pcounter_tc1 + i] = REG_VEC2[i];
								}
								valid_vec2 = false;
								num_sram_access += MULTIPLIER_NUM;
							}

							compute_latency = cycle - prev_cycle; //recore the latency of computation
							if (dram_access_latency > compute_latency)
								cycle = prev_cycle +
										dram_access_latency; //if dram latency costs more cycles, update current cycle using dram latency

							valid_submat_o = true;
							prev_pcounter_n1 = pcounter_n1;
							prev_pcounter_c1 = pcounter_c1;
							prev_pcounter_m1 = pcounter_m1;
							prev_pcounter_n0 = pcounter_n0;
							prev_pcounter_c0 = pcounter_c0;
							prev_pcounter_m0 = pcounter_m0;
						}
					}
				}

				// perform related GCN operation
				/*
				for (pcounter_n1 = 0; pcounter_n1 < DIM_N; pcounter_n1 += TILE_SIZE_TN1) {
					for (pcounter_c1 = 0; pcounter_c1 < DIM_C; pcounter_c1 += TILE_SIZE_TC1) {
						for (pcounter_m1 = 0; pcounter_m1 < DIM_M; pcounter_m1 += TILE_SIZE_TM1) {
							prev_num_dram_access = num_dram_access;

							// load SUBMAT_A
							if (pcounter_m1 != prev_pcounter_m1 || pcounter_n1 != prev_pcounter_n1) {
								
								int num_row = std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1);
								int num_col = std::min(TILE_SIZE_TN1, DIM_C - pcounter_n1);
								generateVirtualData(SUBMAT_A, num_row, num_col, DENSITY_A);

								num_dram_access += countDramAccess(SUBMAT_A);
							}

							// load SUBMAT_BI: uses the power-law pre-loading strategy
							// cout << pcounter_n1 << " " << DIM_M << endl;
							if (pcounter_n1 < 0.2 * DIM_M) {
								// the hot vertices are already pre-loaded in the the on-chip buffer
								if (pcounter_n1 != prev_pcounter_n1 || pcounter_c1 != prev_pcounter_c1) {
									
									int num_row = std::min(TILE_SIZE_TN1, DIM_N - pcounter_n1);
									int num_col = std::min(TILE_SIZE_TC1, DIM_C - pcounter_c1);
									generateVirtualData(SUBMAT_BI, num_row, num_col, 1);
								}
							} else {
								// needed vertices are not pre-loaded in the on-chip buffer
								if (pcounter_n1 != prev_pcounter_n1 || pcounter_c1 != prev_pcounter_c1) {
									int num_row = std::min(TILE_SIZE_TN1, DIM_N - pcounter_n1);
									int num_col = std::min(TILE_SIZE_TC1, DIM_C - pcounter_c1);
									generateVirtualData(SUBMAT_BI, num_row, num_col, 1);

									num_dram_access += countDramAccess(SUBMAT_BI);
								}
							}
							// if (pcounter_n1 != prev_pcounter_n1 || pcounter_c1 != prev_pcounter_c1) {
								
							// 	int num_row = std::min(TILE_SIZE_TN1, DIM_N - pcounter_n1);
							// 	int num_col = std::min(TILE_SIZE_TC1, DIM_C - pcounter_c1);
							// 	generateVirtualData(SUBMAT_BI, num_row, num_col, 1);

							// 	num_dram_access += countDramAccess(SUBMAT_BI);
							// }

							// load SUBMAT_O
							if (pcounter_m1 != prev_pcounter_m1 || pcounter_c1 != prev_pcounter_c1) {
								// check if the SUBMAT_O is valid to write back to DRAM
								if (valid_submat_o && prev_pcounter_m1 != -1 && prev_pcounter_c1 != -1) {
									num_dram_access += countDramAccess(SUBMAT_O);
									valid_submat_o = false;
								}

								int num_row = min(TILE_SIZE_TM1, DIM_M - pcounter_m1);
								int num_col = min(TILE_SIZE_TC1, DIM_C - pcounter_c1);
								generateVirtualData(SUBMAT_O, num_row, num_col, 1);

								num_dram_access += countDramAccess(SUBMAT_O);
							}

							dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access - prev_num_dram_access) / (DRAM_BANDWIDTH / CLOCK_RATE);  // calculate the dram latency

							// on-chip computation
							prev_cycle = cycle;
							
							for (pcounter_tm1 = 0; pcounter_tm1 < std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1); pcounter_tm1++) {
								for (pcounter_tn1 = 0; pcounter_tn1 < std::min(TILE_SIZE_TN1, DIM_N - pcounter_n1); pcounter_tn1++) {
									if (SUBMAT_A[pcounter_tm1][pcounter_tn1] != 0) {
										if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tn1 != prev_pcounter_tn1) {
											REG_SCALAR = SUBMAT_A[pcounter_tm1][pcounter_tn1];
											num_sram_access += 1;
										}

										for (pcounter_tc1 = 0; pcounter_tc1 < std::min(TILE_SIZE_TC1, DIM_C - pcounter_c1); pcounter_tc1 += MULTIPLIER_NUM) {
											// load element of BI to registers
											if (pcounter_tn1 != prev_pcounter_tn1 || pcounter_tc1 != prev_pcounter_tc1) {
												for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC1 - pcounter_tc1); i++) {
													REG_VEC1[i] = SUBMAT_BI[pcounter_tn1][pcounter_tc1 + i];
												}
												num_sram_access += MULTIPLIER_NUM;
											}

											// load elements of O to registers
											if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tc1 != prev_pcounter_tc1) {
												// if replaced, the results in registers need write back
												if (valid_vec2 && prev_pcounter_tm1 != -1 && prev_pcounter_tc1 != -1) {
													for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC1 - prev_pcounter_tc1); i++) {
														SUBMAT_O[prev_pcounter_tm1][prev_pcounter_tc1 + i] = REG_VEC2[i];
													}
													valid_vec2 = false;
													num_sram_access += MULTIPLIER_NUM;
												}

												for (int i = 0;
													i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC1 - pcounter_tc1); i++) {
													REG_VEC2[i] = SUBMAT_O[pcounter_tm1][pcounter_tc1 + i];
												}
												num_sram_access += MULTIPLIER_NUM;
											}
											compute();
											cycle++;
											//printf(" --------------------------------------------current cycle is---%10d-----------------------------------------------\n", cycle);
											valid_vec2 = true;
											prev_pcounter_tc1 = pcounter_tc1;
											prev_pcounter_tn1 = pcounter_tn1;
											prev_pcounter_tm1 = pcounter_tm1;
										}
									}
								}
							}
							if (valid_vec2 && prev_pcounter_tm1 != -1 && prev_pcounter_tc1 != -1) {
								for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC1 - prev_pcounter_tc1); i++) {
									SUBMAT_O[prev_pcounter_tm1][prev_pcounter_tc1 + i] = REG_VEC2[i];
								}
								valid_vec2 = false;
								num_sram_access += MULTIPLIER_NUM;
							}

							compute_latency = cycle - prev_cycle; //recore the latency of computation
							if (dram_access_latency > compute_latency)
								cycle = prev_cycle + dram_access_latency; //if dram latency costs more cycles, update current cycle using dram latency

							valid_submat_o = true;
							prev_pcounter_n1 = pcounter_n1;
							prev_pcounter_c1 = pcounter_c1;
							prev_pcounter_m1 = pcounter_m1;
						}
					}
				}
				*/
				//  store SUBMAT_O
				// if SUBMAT_O is valid, which means it has to write back to DRAM
				if (valid_submat_o && prev_pcounter_m1 != -1 && prev_pcounter_c1 != -1) {
					//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - prev_pcounter_m) * std::min(TILE_SIZE_TC1, DIM_C - prev_pcounter_c1);
					num_dram_access += countDramAccess(SUBMAT_O);
				}
				valid_submat_o = false;

			} else {
				// perform the aggregation first and then combination

				// off-chip communication: B=A*X
				// pre-load data based on the power-law theory
				// calculate the cycle to traverse the CSR to sort the hot vertices
				cycle += DIM_M;

				// for (pcounter_m1 = 0; pcounter_m1 < DIM_M; pcounter_m1 += TILE_SIZE_TM1) {
					// for (pcounter_k1 = 0; pcounter_k1 < DIM_K; pcounter_k1 += TILE_SIZE_TK1) {
						// for (pcounter_n1 = 0; pcounter_n1 < DIM_N; pcounter_n1 < TILE_SIZE_TN1) {


				for (pcounter_n1 = 0; pcounter_n1 < DIM_N; pcounter_n1 += TILE_SIZE_TN1) {
					for (pcounter_k1 = 0; pcounter_k1 < DIM_K; pcounter_k1 += TILE_SIZE_TK1) {
						for (pcounter_m1 = 0; pcounter_m1 < DIM_M; pcounter_m1 < TILE_SIZE_TM1) {			
							prev_num_dram_access = num_dram_access; // record the current number of dram accesses

							// load SUBMAT_A
							if (pcounter_m1 != prev_pcounter_m1 || pcounter_n1 != prev_pcounter_n1) {
								int num_row = std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1);
								int num_col = std::min(TILE_SIZE_TN1, DIM_C - pcounter_n1);
								generateVirtualData(SUBMAT_A, num_row, num_col, DENSITY_A);

								num_dram_access += countDramAccess(SUBMAT_A);
							}

							// load SUBMAT_X
							if (pcounter_n1 < DIM_N) {
								// the needed data are already pre-loaded to the on-chip buffer
								if (pcounter_n1 != prev_pcounter_n1 || pcounter_k1 != prev_pcounter_k1) {
									int num_row = std::min(TILE_SIZE_TN1, DIM_N - pcounter_n1);
									int num_col = std::min(TILE_SIZE_TK1, DIM_K - pcounter_k1);
									generateVirtualData(SUBMAT_A, num_row, num_col, DENSITY_X);
								}
							} else {
								// the needed data is not pre-loaded to the on-chip buffer
								if (pcounter_n1 != prev_pcounter_n1 || pcounter_k1 != prev_pcounter_k1) {
									int num_row = std::min(TILE_SIZE_TN1, DIM_N - pcounter_n1);
									int num_col = std::min(TILE_SIZE_TK1, DIM_K - pcounter_k1);
									generateVirtualData(SUBMAT_A, num_row, num_col, DENSITY_X);

									num_dram_access += countDramAccess(SUBMAT_X);
								}
							}

							// load SUBMAT_BO
							if (pcounter_m1 != prev_pcounter_m1 || pcounter_k1 != prev_pcounter_k1) {
								// if SUBMAT_BO is valid, it has to write back to DRAM
								if (valid_submat_bo && prev_pcounter_m1 != -1 && prev_pcounter_k1 != -1) {
									num_dram_access += countDramAccess(SUBMAT_BO);
									valid_submat_bo = false;
								}

								int num_row = min(TILE_SIZE_TM1, DIM_M - pcounter_m1);
								int num_col = min(TILE_SIZE_TK1, DIM_K - pcounter_k1);
								generateVirtualData(SUBMAT_BO, num_row, num_col, 1);

								num_dram_access += countDramAccess(SUBMAT_BO);
							}

							dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access - prev_num_dram_access) / (DRAM_BANDWIDTH / CLOCK_RATE);

							// on-chip computation
							prev_cycle = cycle;
							for (pcounter_tm1 = 0; pcounter_tm1 < std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1); pcounter_tm1++) {
								for (pcounter_tn1 = 0; pcounter_tn1 < std::min(TILE_SIZE_TN1, DIM_N - pcounter_n1); pcounter_tn1++) {
									if (SUBMAT_A[pcounter_tm1][pcounter_tn1] != 0) {
										if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tn1 != prev_pcounter_tn1) {
											REG_SCALAR = SUBMAT_A[pcounter_tm1][pcounter_tn1];
											num_sram_access += 1;
										}

										for (pcounter_tk1 = 0; pcounter_tk1 < min(TILE_SIZE_TK1, DIM_K - pcounter_k1); pcounter_tk1 += MULTIPLIER_NUM) {
											// load element of X to registers
											if (pcounter_tn1 != prev_pcounter_tn1 || pcounter_tk1 != prev_pcounter_tk1) {
												for (int i = 0; i < min(MULTIPLIER_NUM, TILE_SIZE_TK1 - pcounter_tk1); i++) {
													REG_VEC1[i] = SUBMAT_X[pcounter_tn1][pcounter_tk1 + i];
												}

												num_sram_access += MULTIPLIER_NUM;
											}

											// load element of B to registers
											if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tk1 != prev_pcounter_tk1) {
												if (valid_vec2 && prev_pcounter_tm1 != -1 && prev_pcounter_tk1 != -1) {
													for (int i = 0; i < min(MULTIPLIER_NUM, TILE_SIZE_TK1 - prev_pcounter_tk1); i++) {
														SUBMAT_BO[prev_pcounter_tm1][prev_pcounter_tk1 + i] = REG_VEC2[i];
													}
													valid_vec2 = false;
													num_sram_access += MULTIPLIER_NUM;
												}

												for (int i = 0; i < min(MULTIPLIER_NUM, TILE_SIZE_TK1 - pcounter_tk1); i++) {
													REG_VEC2[i] = SUBMAT_BO[pcounter_tm1][pcounter_tk1 + i];
												}

												num_sram_access += MULTIPLIER_NUM;
											}

											compute();
											cycle = cycle + 0.9;
											valid_vec2 = true;
											prev_pcounter_tm1 = pcounter_tm1;
											prev_pcounter_tk1 = pcounter_tk1;
											prev_pcounter_tn1 = pcounter_tn1;
										}
									}
								}
							}

							if (valid_vec2 && prev_pcounter_tm1 != -1 && prev_pcounter_tk1 != -1) {
								for (int i = 0; i < min(MULTIPLIER_NUM, TILE_SIZE_TK1 -prev_pcounter_tk1); i++) {
									SUBMAT_BO[prev_pcounter_tm1][pcounter_tk1 + i] = REG_VEC2[i];
								}
								valid_vec2 = false;
								num_sram_access += MULTIPLIER_NUM;
							}

							compute_latency = cycle - prev_cycle;

							if (dram_access_latency > compute_latency) {
								cycle = prev_cycle + dram_access_latency;
							}

							valid_submat_bo = true;
							prev_pcounter_m1 = pcounter_m1;
							prev_pcounter_n1 = pcounter_n1;
							prev_pcounter_k1 = pcounter_k1;
						}
					}
				}

				if (valid_submat_bo && prev_pcounter_m1 != -1 && prev_pcounter_k1 != -1) {
					num_dram_access += countDramAccess(SUBMAT_BO);
				}

				valid_submat_bo = false;
				
				// off-chip communication: O=B*W
				for (pcounter_m0 = 0; pcounter_m0 < DIM_M; pcounter_m0 += TILE_SIZE_TM0) {
					for (pcounter_c0 = 0; pcounter_c0 < DIM_C; pcounter_c0 += TILE_SIZE_TC0) {
						for (pcounter_k0 = 0; pcounter_k0 < DIM_K; pcounter_k0 += TILE_SIZE_TK0) {
							prev_num_dram_access = num_dram_access;

							// load SUBMAT_BI
							if (pcounter_m0 != prev_pcounter_m0 || pcounter_k0 != prev_pcounter_k0) {
								int num_row = min(TILE_SIZE_TM0, DIM_M - pcounter_m0);
								int num_col = min(TILE_SIZE_TK0, DIM_K - pcounter_k0);
								generateVirtualData(SUBMAT_BI, num_row, num_col, 1);

								num_dram_access += countDramAccess(SUBMAT_BI);
							}

							// load SUBMAT_W
							if (pcounter_k0 != prev_pcounter_k0 || pcounter_c0 != prev_pcounter_c0) {
								int num_row = min(TILE_SIZE_TK0, DIM_K - pcounter_k0);
								int num_col = min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
								generateVirtualData(SUBMAT_W, num_row, num_col, DENSITY_W);

								num_dram_access += countDramAccess(SUBMAT_W);
							}

							// load SUBMAT_O
							if (pcounter_m0 != prev_pcounter_m0 || pcounter_c0 != prev_pcounter_c0) {
								// if replace or not
								if (valid_submat_o && prev_pcounter_m0 != -1 && prev_pcounter_c0 != -1) {
									num_dram_access += countDramAccess(SUBMAT_O);
									valid_submat_o = false;
								}

								int num_row = min(TILE_SIZE_TM0, DIM_M - pcounter_m0);
								int num_col = min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
								generateVirtualData(SUBMAT_O, num_row, num_col, 1);

								num_dram_access += countDramAccess(SUBMAT_O);
							}

							dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access - prev_num_dram_access) / (DRAM_BANDWIDTH / CLOCK_RATE);

							// on-chip computation
							prev_cycle = cycle;
							for (pcounter_tm0 = 0; pcounter_tm0 < min(TILE_SIZE_TM0, DIM_M - pcounter_m0); pcounter_tm0++) {
								for (pcounter_tk0 = 0; pcounter_tk0 < min(TILE_SIZE_TK0, DIM_K - pcounter_k0); pcounter_tk0++) {
									if (SUBMAT_BI[pcounter_tm0][pcounter_tk0] != 0) {
										if (pcounter_tm0 != prev_pcounter_tm0 || pcounter_tk0 != prev_pcounter_tk0) {
											REG_SCALAR = SUBMAT_BI[pcounter_tm0][pcounter_tk0];
											num_sram_access += 1;
										}

										for (pcounter_tc0 = 0; pcounter_tc0 < min(TILE_SIZE_TC0, DIM_C - pcounter_c0); pcounter_tc0 += MULTIPLIER_NUM) {
											// load W to registers
											if (pcounter_tk0 != prev_pcounter_tk0 || pcounter_tc0 != prev_pcounter_tc0) {
												for (int i = 0; i < min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc0); i++) {
													REG_VEC1[i] = SUBMAT_W[pcounter_tk0][pcounter_tc0 + i];
												}

												num_sram_access += MULTIPLIER_NUM;
											}

											// load O to registers
											if (pcounter_tm0 != prev_pcounter_tm0 || pcounter_tc0 != prev_pcounter_tc0) {
												// check if replace
												if (valid_vec1 && prev_pcounter_tm0 != -1 && prev_pcounter_tc0 != 01) {
													for (int i = 0; i < min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc0); i++) {
														SUBMAT_O[pcounter_tm0][pcounter_tc0 + i] = REG_VEC2[i];
													}
													
													valid_vec1 = false;
													num_sram_access += MULTIPLIER_NUM;
												}

												for (int i = 0; i < min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc0); i++) {
													REG_VEC2[i] = SUBMAT_O[pcounter_tm0][pcounter_tc0 + i];
												}

												num_sram_access += MULTIPLIER_NUM;
											}

											compute();
											cycle += 0.9;
											valid_vec1 = true;

											prev_pcounter_tm0 = pcounter_tm0;
											prev_pcounter_tk0 = pcounter_tk0;
											prev_pcounter_tc0 = pcounter_tc0;
										}
									}
								}
							}

							// store SUBMAT_O
							if (valid_vec1 && prev_pcounter_tm0 != -1 && prev_pcounter_tc0 != -1) {
								for (int i = 0; i < min(MULTIPLIER_NUM, TILE_SIZE_TC0 - prev_pcounter_tc0); i++) {
									SUBMAT_O[pcounter_tm0][pcounter_tc0 + i] = REG_VEC2[i];
								}
								valid_vec1 = false;
								num_sram_access += MULTIPLIER_NUM;
							}

							compute_latency = cycle - prev_cycle;
							if (dram_access_latency > compute_latency) {
								cycle = prev_cycle + dram_access_latency;
							}
							
							valid_submat_o = true;
							prev_pcounter_c0 = pcounter_c0;
							prev_pcounter_k0 = pcounter_k0;
							prev_pcounter_m0 = pcounter_m0;
						}
					}
				}

				if (valid_submat_o && prev_pcounter_m0 != -1 && prev_pcounter_c0 != -1) {
					num_dram_access += countDramAccess(SUBMAT_O);
				}
				valid_submat_o = false;
			}
		} else {
			std::cout << "Wrong Accelerator Type!" << std::endl;
		}
	} else {
		if (!LOOP_FUSION) {
			prev_pcounter_n0 = -1;
			prev_pcounter_c0 = -1;
			prev_pcounter_k0 = -1;
			prev_pcounter_tn0 = -1;
			prev_pcounter_tc0 = -1;
			prev_pcounter_tk0 = -1;
			prev_pcounter_m1 = -1;
			prev_pcounter_c1 = -1;
			prev_pcounter_n1 = -1;
			prev_pcounter_tm1 = -1;
			prev_pcounter_tc1 = -1;
			prev_pcounter_tn1 = -1;

			//off-chip communication: SPMM1, B=X*W
			//load matrix X, W, B to SUBMAT_X, SUBMAT_W, SUBMAT_BO

			for (pcounter_n0 = 0; pcounter_n0 < DIM_N; pcounter_n0 += TILE_SIZE_TN0) {
				for (pcounter_c0 = 0; pcounter_c0 < DIM_C; pcounter_c0 += TILE_SIZE_TC0) {
					for (pcounter_k0 = 0; pcounter_k0 < DIM_K; pcounter_k0 += TILE_SIZE_TK0) {
						prev_num_dram_access = num_dram_access; //record the current number of dram accesses

						if (pcounter_n0 != prev_pcounter_n0 || pcounter_k0 != prev_pcounter_k0) {
							//  load SUBMAT_X
							int num_row = std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
							int num_col = std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0);
							generateVirtualData(SUBMAT_X, num_row, num_col, DENSITY_X);

							//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TK, DIM_K - pcounter_k);
							num_dram_access += countDramAccess(SUBMAT_X);
						}

						if (pcounter_k0 != prev_pcounter_k0 || pcounter_c0 != prev_pcounter_c0) {
							//  load SUBMAT_W
							int num_row = std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0);
							int num_col = std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							generateVirtualData(SUBMAT_W, num_row, num_col, DENSITY_W);

							//num_dram_access += std::min(TILE_SIZE_TK, DIM_K - pcounter_k) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							num_dram_access += countDramAccess(SUBMAT_W);
						}

						if (pcounter_n0 != prev_pcounter_n0 || pcounter_c0 != prev_pcounter_c0) {
							//  load SUBMAT_BO
							// if SUBMAT_BO is valid, which means it has to write back to DRAM
							if (valid_submat_bo && prev_pcounter_n0 != -1 && prev_pcounter_c0 != -1) {
								//store SUBMAT_BO to DRAM

								//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
								num_dram_access += countDramAccess(SUBMAT_BO);
								valid_submat_bo = false;
							}

							//  load SUBMAT_BO
							int num_row = std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0);
							int num_col = std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							generateVirtualData(SUBMAT_BO, num_row, num_col, 1);

							//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							num_dram_access += countDramAccess(SUBMAT_BO);
						}

						dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access - prev_num_dram_access) /
																   (DRAM_BANDWIDTH /
																	CLOCK_RATE); //calculate the dram access latency

						//on-chip computation
						// loop order tn0->tk->tc0
						prev_cycle = cycle;
						for (pcounter_tn0 = 0;
							 pcounter_tn0 < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); pcounter_tn0++) {
							for (pcounter_tk0 = 0;
								 pcounter_tk0 < std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0); pcounter_tk0++) {
								if (SUBMAT_X[pcounter_tn0][pcounter_tk0] != 0) {
									if (pcounter_tn0 != prev_pcounter_tn0 || pcounter_tk0 != prev_pcounter_tk0) {
										REG_SCALAR = SUBMAT_X[pcounter_tn0][pcounter_tk0];
										num_sram_access += 1;
									}

									for (pcounter_tc0 = 0; pcounter_tc0 < std::min(TILE_SIZE_TC0, DIM_C -
																								  pcounter_c0); pcounter_tc0 += MULTIPLIER_NUM) {
										// load element of W to registers
										if (pcounter_tk0 != prev_pcounter_tk0 || pcounter_tc0 != prev_pcounter_tc0) {
											for (int i = 0;
												 i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc0); i++) {
												REG_VEC1[i] = SUBMAT_W[pcounter_tk0][pcounter_tc0 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}

										// load elements of B to registers
										if (pcounter_tn0 != prev_pcounter_tn0 || pcounter_tc0 != prev_pcounter_tc0) {
											// if replaced, the results in registers need write back
											if (valid_vec2 && prev_pcounter_tn0 != -1 && prev_pcounter_tc0 != -1) {
												for (int i = 0; i < std::min(MULTIPLIER_NUM,
																			 TILE_SIZE_TC0 - prev_pcounter_tc0); i++) {
													SUBMAT_BO[prev_pcounter_tn0][prev_pcounter_tc0 + i] = REG_VEC2[i];
												}
												valid_vec2 = false;
												num_sram_access += MULTIPLIER_NUM;
											}

											for (int i = 0;
												 i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc0); i++) {
												REG_VEC2[i] = SUBMAT_BO[pcounter_tn0][pcounter_tc0 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}
										compute();
										cycle++;
										//printf(" --------------------------------------------current cycle is---%10d-----------------------------------------------\n", cycle);
										valid_vec2 = true;
										prev_pcounter_tc0 = pcounter_tc0;
										prev_pcounter_tk0 = pcounter_tk0;
										prev_pcounter_tn0 = pcounter_tn0;
									}
								}
							}
						}
						if (valid_vec2 && prev_pcounter_tn0 != -1 && prev_pcounter_tc0 != -1) {
							for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - prev_pcounter_tc0); i++) {
								SUBMAT_BO[prev_pcounter_tn0][prev_pcounter_tc0 + i] = REG_VEC2[i];
							}
							valid_vec2 = false;
							num_sram_access += MULTIPLIER_NUM;
						}

						compute_latency = cycle - prev_cycle; //recore the latency of computation
						if (dram_access_latency > compute_latency)
							cycle = prev_cycle +
									dram_access_latency; //if dram latency costs more cycles, update current cycle using dram latency

						prev_pcounter_k0 = pcounter_k0;
						valid_submat_bo = true;
						prev_pcounter_c0 = pcounter_c0;
						prev_pcounter_n0 = pcounter_n0;
					}
				}
			}
			//  store SUBMAT_BO
			// if SUBMAT_BO is valid, which means it has to write back to DRAM
			if (valid_submat_bo && prev_pcounter_n0 != -1 && prev_pcounter_c0 != -1) {
				//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - prev_pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - prev_pcounter_c0);
				num_dram_access += countDramAccess(SUBMAT_BO);
			}
			valid_submat_bo = false;

			//off-chip communication: SPMM2, O=A*B
			//load matrix A, B, O to SUBMAT_A, SUBMAT_B, SUBMAT_O
			for (pcounter_m1 = 0; pcounter_m1 < DIM_M; pcounter_m1 += TILE_SIZE_TM1) {
				for (pcounter_c1 = 0; pcounter_c1 < DIM_C; pcounter_c1 += TILE_SIZE_TC1) {
					for (pcounter_n1 = 0; pcounter_n1 < DIM_N; pcounter_n1 += TILE_SIZE_TN1) {
						prev_num_dram_access = num_dram_access; //record the current number of dram accesses

						if (pcounter_m1 != prev_pcounter_m1 || pcounter_n1 != prev_pcounter_n1) {
							//  load SUBMAT_A
							int num_row = std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1);
							int num_col = std::min(TILE_SIZE_TN1, DIM_C - pcounter_n1);
							generateVirtualData(SUBMAT_A, num_row, num_col, DENSITY_A);

							//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TN1, DIM_N - pcounter_n1);
							num_dram_access += countDramAccess(SUBMAT_A);
						}
						
						if (pcounter_n1 != prev_pcounter_n1 || pcounter_c1 != prev_pcounter_c1) {
							//  load SUBMAT_BI
							int num_row = std::min(TILE_SIZE_TN1, DIM_N - pcounter_n1);
							int num_col = std::min(TILE_SIZE_TC1, DIM_C - pcounter_c1);
							generateVirtualData(SUBMAT_BI, num_row, num_col, 1);

							//num_dram_access += std::min(TILE_SIZE_TN1, DIM_N - pcounter_n1) * std::min(TILE_SIZE_TC1, DIM_C - pcounter_c1);
							num_dram_access += countDramAccess(SUBMAT_BI);
						}

						if (pcounter_m1 != prev_pcounter_m1 || pcounter_c1 != prev_pcounter_c1) {
							//  store SUBMAT_O
							// if SUBMAT_BO is valid, which means it has to write back to DRAM
							if (valid_submat_o && prev_pcounter_m1 != -1 && prev_pcounter_c1 != -1) {
								//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TC1, DIM_C - pcounter_c1);
								num_dram_access += countDramAccess(SUBMAT_O);
								valid_submat_o = false;
							}

							//  load SUBMAT_O
							int num_row = std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1);
							int num_col = std::min(TILE_SIZE_TC1, DIM_C - pcounter_c1);
							generateVirtualData(SUBMAT_O, num_row, num_col, 1);

							//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TC1, DIM_C - pcounter_c1);
							num_dram_access += countDramAccess(SUBMAT_O);
						}

						//cout << num_dram_access << endl;
						dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access - prev_num_dram_access) /
																   (DRAM_BANDWIDTH /
																	CLOCK_RATE); //calculate the dram access latency

						//on-chip computation
						// loop order tn0->tk->tc0
						prev_cycle = cycle;
						for (pcounter_tm1 = 0;
							 pcounter_tm1 < std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1); pcounter_tm1++) {
							for (pcounter_tn1 = 0;
								 pcounter_tn1 < std::min(TILE_SIZE_TN1, DIM_N - pcounter_n1); pcounter_tn1++) {
								if (SUBMAT_A[pcounter_tm1][pcounter_tn1] != 0) {
									if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tn1 != prev_pcounter_tn1) {
										REG_SCALAR = SUBMAT_A[pcounter_tm1][pcounter_tn1];
										num_sram_access += 1;
									}

									for (pcounter_tc1 = 0; pcounter_tc1 < std::min(TILE_SIZE_TC1, DIM_C -
																								  pcounter_c1); pcounter_tc1 += MULTIPLIER_NUM) {
										// load element of BI to registers
										if (pcounter_tn1 != prev_pcounter_tn1 || pcounter_tc1 != prev_pcounter_tc1) {
											for (int i = 0;
												 i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC1 - pcounter_tc1); i++) {
												REG_VEC1[i] = SUBMAT_BI[pcounter_tn1][pcounter_tc1 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}

										// load elements of O to registers
										if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tc1 != prev_pcounter_tc1) {
											// if replaced, the results in registers need write back
											if (valid_vec2 && prev_pcounter_tm1 != -1 && prev_pcounter_tc1 != -1) {
												for (int i = 0; i < std::min(MULTIPLIER_NUM,
																			 TILE_SIZE_TC1 - prev_pcounter_tc1); i++) {
													SUBMAT_O[prev_pcounter_tm1][prev_pcounter_tc1 + i] = REG_VEC2[i];
												}
												valid_vec2 = false;
												num_sram_access += MULTIPLIER_NUM;
											}

											for (int i = 0;
												 i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC1 - pcounter_tc1); i++) {
												REG_VEC2[i] = SUBMAT_O[pcounter_tm1][pcounter_tc1 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}
										compute();
										cycle++;
										//printf(" --------------------------------------------current cycle is---%10d-----------------------------------------------\n", cycle);
										valid_vec2 = true;
										prev_pcounter_tc1 = pcounter_tc1;
										prev_pcounter_tn1 = pcounter_tn1;
										prev_pcounter_tm1 = pcounter_tm1;
									}
								}
							}
						}
						if (valid_vec2 && prev_pcounter_tm1 != -1 && prev_pcounter_tc1 != -1) {
							for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC1 - prev_pcounter_tc1); i++) {
								SUBMAT_O[prev_pcounter_tm1][prev_pcounter_tc1 + i] = REG_VEC2[i];
							}
							valid_vec2 = false;
							num_sram_access += MULTIPLIER_NUM;
						}

						compute_latency = cycle - prev_cycle; //recore the latency of computation
						if (dram_access_latency > compute_latency)
							cycle = prev_cycle +
									dram_access_latency; //if dram latency costs more cycles, update current cycle using dram latency

						valid_submat_o = true;
						prev_pcounter_n1 = pcounter_n1;
						prev_pcounter_c1 = pcounter_c1;
						prev_pcounter_m1 = pcounter_m1;
					}
				}
			}

			//  store SUBMAT_O
			// if SUBMAT_O is valid, which means it has to write back to DRAM
			if (valid_submat_o && prev_pcounter_m1 != -1 && prev_pcounter_c1 != -1) {
				//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - prev_pcounter_m) * std::min(TILE_SIZE_TC1, DIM_C - prev_pcounter_c1);
				num_dram_access += countDramAccess(SUBMAT_O);
			}
			valid_submat_o = false;
		}
			// Enable loop fusion
		else {
			prev_pcounter_n0 = -1;
			prev_pcounter_c0 = -1;
			prev_pcounter_k0 = -1;
			prev_pcounter_tn0 = -1;
			prev_pcounter_tc0 = -1;
			prev_pcounter_tk0 = -1;

			prev_pcounter_m1 = -1;
			prev_pcounter_c1 = -1;
			prev_pcounter_n1 = -1;
			prev_pcounter_tm1 = -1;
			prev_pcounter_tc1 = -1;
			prev_pcounter_tn1 = -1;

			for (pcounter_n0 = 0; pcounter_n0 < DIM_N; pcounter_n0 += TILE_SIZE_TN0) {
				for (pcounter_c0 = 0; pcounter_c0 < DIM_C; pcounter_c0 += TILE_SIZE_TC0) {
					for (pcounter_k0 = 0; pcounter_k0 < DIM_K; pcounter_k0 += TILE_SIZE_TK0) {
						prev_num_dram_access = num_dram_access; //record the current number of dram accesses

						if (pcounter_n0 != prev_pcounter_n0 || pcounter_k0 != prev_pcounter_k0) {
							//  load SUBMAT_X
							int num_row = std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
							int num_col = std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0);
							generateVirtualData(SUBMAT_X, num_row, num_col, DENSITY_X);

							num_dram_access += countDramAccess(SUBMAT_X);
						}

						if (pcounter_k0 != prev_pcounter_k0 || pcounter_c0 != prev_pcounter_c0) {
							//  load SUBMAT_W
							int num_row = std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0);
							int num_col = std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							generateVirtualData(SUBMAT_W, num_row, num_col, DENSITY_W);

							//num_dram_access += std::min(TILE_SIZE_TK, DIM_K - pcounter_k) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							num_dram_access += countDramAccess(SUBMAT_W);
						}

						if (pcounter_n0 != prev_pcounter_n0 || pcounter_c0 != prev_pcounter_c0) {
							//  load SUBMAT_BO
							// if SUBMAT_BO is valid, which means it has to write back to DRAM
							if (valid_submat_bo && prev_pcounter_n0 != -1 && prev_pcounter_c0 != -1) {
								//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
								//num_dram_access += countDramAccess(SUBMAT_X);
								valid_submat_bo = false;
							}

							//  load SUBMAT_BO
							int num_row = std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
							int num_col = std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							generateVirtualData(SUBMAT_BO, num_row, num_col, 1);

							//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
						}

						dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access - prev_num_dram_access) /
																   (DRAM_BANDWIDTH /
																	CLOCK_RATE); //calculate the dram access latency

						//on-chip computation
						// loop order tn0->tk->tc0
						for (pcounter_tn0 = 0;
							 pcounter_tn0 < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); pcounter_tn0++) {
							for (pcounter_tk0 = 0;
								 pcounter_tk0 < std::min(TILE_SIZE_TK0, DIM_K - pcounter_k0); pcounter_tk0++) {
								if (SUBMAT_X[pcounter_tn0][pcounter_tk0] != 0) {
									if (pcounter_tn0 != prev_pcounter_tn0 || pcounter_tk0 != prev_pcounter_tk0) {
										REG_SCALAR = SUBMAT_X[pcounter_tn0][pcounter_tk0];
										num_sram_access += 1;
									}

									for (pcounter_tc0 = 0; pcounter_tc0 < std::min(TILE_SIZE_TC0, DIM_C -
																								  pcounter_c0); pcounter_tc0 += MULTIPLIER_NUM) {
										// load element of W to registers
										if (pcounter_tk0 != prev_pcounter_tk0 || pcounter_tc0 != prev_pcounter_tc0) {
											for (int i = 0;
												 i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc0); i++) {
												REG_VEC1[i] = SUBMAT_W[pcounter_tk0][pcounter_tc0 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}

										// load elements of B to registers
										if (pcounter_tn0 != prev_pcounter_tn0 || pcounter_tc0 != prev_pcounter_tc0) {
											// if replaced, the results in registers need write back
											if (valid_vec2 && prev_pcounter_tn0 != -1 && prev_pcounter_tc0 != -1) {
												for (int i = 0; i < std::min(MULTIPLIER_NUM,
																			 TILE_SIZE_TC0 - prev_pcounter_tc0); i++) {
													SUBMAT_BO[prev_pcounter_tn0][prev_pcounter_tc0 + i] = REG_VEC2[i];
												}
												valid_vec2 = false;
												num_sram_access += MULTIPLIER_NUM;
											}

											for (int i = 0;
												 i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc0); i++) {
												REG_VEC2[i] = SUBMAT_BO[pcounter_tn0][pcounter_tc0 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}
										compute();
										cycle++;

										//printf(" --------------------------------------------current cycle is---%10d-----------------------------------------------\n", cycle);
										valid_vec2 = true;
										prev_pcounter_tc0 = pcounter_tc0;
										prev_pcounter_tk0 = pcounter_tk0;
										prev_pcounter_tn0 = pcounter_tn0;
									}
								}
							}
						}
						if (valid_vec2 && prev_pcounter_tn0 != -1 && prev_pcounter_tc0 != -1) {
							for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - prev_pcounter_tc0); i++) {
								SUBMAT_BO[prev_pcounter_tn0][prev_pcounter_tc0 + i] = REG_VEC2[i];
							}
							valid_vec2 = false;
							num_sram_access += MULTIPLIER_NUM;
						}

						compute_latency = cycle - prev_cycle; //recore the latency of computation
						if (dram_access_latency > compute_latency)
							cycle = prev_cycle +
									dram_access_latency; //if dram latency costs more cycles, update current cycle using dram latency

						prev_pcounter_k0 = pcounter_k0;
						valid_submat_bo = true;
						prev_pcounter_c0 = pcounter_c0;
						prev_pcounter_n0 = pcounter_n0;
					}

					cout << num_dram_access << endl;
					for (pcounter_m1 = 0; pcounter_m1 < DIM_M; pcounter_m1 += TILE_SIZE_TM1) {
						prev_num_dram_access = num_dram_access; //record the current number of dram accesses

						// off-chip communication, load SUBMAT_A, SUBMAT_O
						if (pcounter_m1 != prev_pcounter_m1 || pcounter_n0 != prev_pcounter_n0) {
							//  load SUBMAT_A
							int num_row = std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1);
							int num_col = std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
							generateVirtualData(SUBMAT_A, num_row, num_col, DENSITY_A);

							//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
							num_dram_access += countDramAccess(SUBMAT_A);
						}

						if (pcounter_n0 != prev_pcounter_n0 || pcounter_c0 != prev_pcounter_c0) {
							//  load SUBMAT_BI
							int num_row = std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0);
							int num_col = std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							generateVirtualData(SUBMAT_BI, num_row, num_col, 1);

							//num_dram_access += std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							num_dram_access += countDramAccess(SUBMAT_BI);
						}

						if (pcounter_m1 != prev_pcounter_m1 || pcounter_c0 != prev_pcounter_c0) {
							//  store SUBMAT_O
							// if SUBMAT_BO is valid, which means it has to write back to DRAM
							if (valid_submat_o && prev_pcounter_m1 != -1 && prev_pcounter_c0 != -1) {
								//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
								num_dram_access += countDramAccess(SUBMAT_O);
								valid_submat_o = false;
							}

							//  load SUBMAT_O
							int num_row = std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1);
							int num_col = std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							generateVirtualData(SUBMAT_O, num_row, num_col, 1);

							//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - pcounter_m) * std::min(TILE_SIZE_TC0, DIM_C - pcounter_c0);
							num_dram_access += countDramAccess(SUBMAT_O);
						}

						dram_access_latency = DRAM_SETUP_LATENCY + (num_dram_access - prev_num_dram_access) /
																   (DRAM_BANDWIDTH /
																	CLOCK_RATE); //calculate the dram access latency

						//on-chip computation
						// loop order tn0->tk->tc0
						for (pcounter_tm1 = 0;
							 pcounter_tm1 < std::min(TILE_SIZE_TM1, DIM_M - pcounter_m1); pcounter_tm1++) {
							for (pcounter_tn1 = 0;
								 pcounter_tn1 < std::min(TILE_SIZE_TN0, DIM_N - pcounter_n0); pcounter_tn1++) {
								if (SUBMAT_A[pcounter_tm1][pcounter_tn1] != 0) {
									if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tn1 != prev_pcounter_tn1) {
										REG_SCALAR = SUBMAT_A[pcounter_tm1][pcounter_tn1];
										num_sram_access += 1;
									}

									for (pcounter_tc1 = 0; pcounter_tc1 < std::min(TILE_SIZE_TC0, DIM_C -
																								  pcounter_c0); pcounter_tc1 += MULTIPLIER_NUM) {
										// load element of BI to registers
										if (pcounter_tn1 != prev_pcounter_tn1 || pcounter_tc1 != prev_pcounter_tc1) {
											for (int i = 0;
												 i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc1); i++) {
												REG_VEC1[i] = SUBMAT_BO[pcounter_tn1][pcounter_tc1 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}

										// load elements of O to registers
										if (pcounter_tm1 != prev_pcounter_tm1 || pcounter_tc1 != prev_pcounter_tc1) {
											// if replaced, the results in registers need write back
											if (valid_vec2 && prev_pcounter_tm1 != -1 && prev_pcounter_tc1 != -1) {
												for (int i = 0; i < std::min(MULTIPLIER_NUM,
																			 TILE_SIZE_TC0 - prev_pcounter_tc1); i++) {
													SUBMAT_O[prev_pcounter_tm1][prev_pcounter_tc1 + i] = REG_VEC2[i];
												}
												valid_vec2 = false;
												num_sram_access += MULTIPLIER_NUM;
											}

											for (int i = 0;
												 i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - pcounter_tc1); i++) {
												REG_VEC2[i] = SUBMAT_O[pcounter_tm1][pcounter_tc1 + i];
											}
											num_sram_access += MULTIPLIER_NUM;
										}
										compute();
										cycle++;
										//printf(" --------------------------------------------current cycle is---%10d-----------------------------------------------\n", cycle);
										valid_vec2 = true;
										prev_pcounter_tc1 = pcounter_tc1;
										prev_pcounter_tn1 = pcounter_tn1;
										prev_pcounter_tm1 = pcounter_tm1;
									}
								}
							}
						}
						if (valid_vec2 && prev_pcounter_tm1 != -1 && prev_pcounter_tc1 != -1) {
							for (int i = 0; i < std::min(MULTIPLIER_NUM, TILE_SIZE_TC0 - prev_pcounter_tc1); i++) {
								SUBMAT_O[prev_pcounter_tm1][prev_pcounter_tc1 + i] = REG_VEC2[i];
							}
							valid_vec2 = false;
							num_sram_access += MULTIPLIER_NUM;
						}

						compute_latency = cycle - prev_cycle; //recore the latency of computation
						if (dram_access_latency > compute_latency)
							cycle = prev_cycle +
									dram_access_latency; //if dram latency costs more cycles, update current cycle using dram latency

						valid_submat_o = true;
						prev_pcounter_n1 = pcounter_n1;
						prev_pcounter_c1 = pcounter_c1;
						prev_pcounter_m1 = pcounter_m1;
						prev_pcounter_n0 = pcounter_n0;
						prev_pcounter_c0 = pcounter_c0;
						prev_pcounter_m0 = pcounter_m0;
					}
				}
			}
			//  store SUBMAT_O
			// if SUBMAT_O is valid, which means it has to write back to DRAM
			if (valid_submat_o && prev_pcounter_m1 != -1 && prev_pcounter_c0 != -1) {
				//num_dram_access += std::min(TILE_SIZE_TM, DIM_M - prev_pcounter_m) * std::min(TILE_SIZE_TC0, DIM_C - prev_pcounter_c1);
				num_dram_access += countDramAccess(SUBMAT_O);
			}
			valid_submat_o = false;
		}
	}

	// reduce cycle because of systolic array

	end_sim = clock();
	double cost_time = (double)(end_sim - begin_sim);

	//save the results to result.txt
	//	writeResultToDram(DATASET);

	printf("\n**********************************************************************************************************************\n");
	printf(" ------------------------------------------------Sim-2 has done-------------------------------------------------------\n");
	printf(" ---------------------------------------Total cycle number is---%10d-----------------------------------------------\n", cycle - 1);
	printf(" --------------------------------------Number of DRAM accesses---%10d---------------------------------------------\n", num_dram_access);
	printf(" --------------------------------------Number of SRAM accesses---%10d---------------------------------------------\n", num_sram_access);
	//FILE *freport_test = fopen("test_report.txt","w+");
	//fprintf(freport_test, "\n\nexecution all time is : %lf seconds.\n\n", cost_time / cycle_TCK);
	//fclose(freport_test);
	printf("**********************************************************************************************************************\n");
}