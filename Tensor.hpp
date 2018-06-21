/*A generic class that defines a tensor for 3 dimensions. A task to-do is to extend this function to a generic dimension. 
 * This class has been written keeping in minded the strain tensor. But this class be used for any material property that is expressed as a tensor
 */ 


class tensor
{
	
	public:
	
	double matrix[3][3];
	tensor()
	{
		
		for(int i = 0 ; i < 3 ; i++)
		{
			for(int j = 0 ; j < 3 ; j++)
			{
				matrix[i][j] = 0.0 ;
			}
		}
	}
};



