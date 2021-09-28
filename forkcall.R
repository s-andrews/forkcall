### Setup and Preferences ####

forward_count_file <- "ExampleData/example_fwd.txt"
reverse_count_file <- "ExampleData/example_rev.txt"

output_dir <- "ExampleData"

min_count <- 25
max_count <- 200

# This value is specified in numbers of observations
# to smooth over
loess_smoothing_level <- 50
min_block_size <- 10

# Draw diagnostic plots to help optimise the setup?
draw_diagnostic_plots <- TRUE


#############################################################

library(tidyverse)

# Create the output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE)

# Define the function to call zones
call_zones <- function(positions,direction_values, sample, chr) {
  
  # Put the data into a tibble
  tibble(
    Position = positions,
    Directionality=direction_values
  ) %>%
    arrange(Position)-> test_data
  
  
  # We need to calculate the proportion of data to
  # smooth over based on the smoothing level config
  # value and the amount of data
  effective_smoothing_level <- loess_smoothing_level/length(positions)
  
  # Smooth with loess smoothing
  test_data %>%
    mutate(loess=predict(loess(Directionality~Position, span=effective_smoothing_level))) -> test_data
  
  # Calculate the first derivative of the smoothed directionality
  test_data %>%
    mutate(derivative=c(0,diff(loess))) -> calling_data
  
  # Call the direction of each transition
  calling_data %>%
    mutate(direction=if_else(derivative>=0, 1, -1)) -> calling_data
  
  # Use RLE to find blocks with the same direction
  rle(calling_data$direction) -> rle_data
  tibble(position_block_size=rle_data[[1]], direction=rle_data[[2]]) -> rle_data
  
  # Annotated the index and true position
  rle_data %>%
    mutate(
      position_index = cumsum(position_block_size),
      position=calling_data$Position[position_index]
    ) -> rle_data
  
  
  # Filter to remove small blocks
  rle_data %>%
    filter(position_block_size >= min_block_size) -> rle_data
  
  # Remove redundant blocks after filtering
  rle_data %>%
    filter(direction != lead(direction, default=TRUE)) -> rle_data
  
  # We need to recalculate the block sizes from the filtered data
  rle_data %>%
    mutate(position_block_size = position_index - lag(position_index,default = 1)) -> rle_data
  
  
  if (draw_diagnostic_plots) {
    rle_data %>%
      filter(position < 25000000) -> rle_for_plot
    
    # Draw a diagnostic plot
    test_data %>%
      filter(Position < 25000000) %>%
      pivot_longer(cols=-Position, names_to="DataType",values_to="value") %>%
      ggplot(aes(x=Position/1000000,y=value,colour=DataType)) +
      geom_line(size=1) +
      theme_bw() +
      geom_vline(xintercept=rle_for_plot$position/1000000,colour="green4")+
      scale_colour_manual(values=c("#555555","red2")) +
      xlab("Position (Mbp)")+
      ggtitle(paste0(sample," chr",chr," smoothing=",loess_smoothing_level," filtering=",min_block_size)) -> diagplot
    
    # Save the plot
    dir.create(paste0(output_dir,"/DiagPlots"),showWarnings = FALSE)
    ggsave(paste0(output_dir,"/DiagPlots/sample",sample,"_chr",chr,".png"),plot=diagplot, width = 14, height=2)
  }  

  # Add additional annotation
  rle_data %>%
    rename(
      end_position=position, 
      end_index=position_index
    ) %>%
    mutate(
      start_position=lag(end_position, default = calling_data$Position[1]),
      start_index=lag(end_index, default = 1),
      value_change=calling_data$loess[end_index]-calling_data$loess[start_index],
      start_value=calling_data$loess[start_index],
      end_value=calling_data$loess[end_index],
      crosses_zero=(start_value>0 & end_value<0) | (start_value<0 & end_value>0),
      chr_distance=end_position-start_position,
      measure_density = (end_index-start_index)/((end_position-start_position)/10000)
    ) %>%
    select(
      direction,
      start_index, 
      end_index, 
      position_block_size,
      start_position, 
      end_position, 
      chr_distance, 
      start_value,
      end_value, 
      value_change, 
      measure_density,
      crosses_zero
    ) -> rle_data
  
  return(rle_data)
  
}


# Read the raw counts
read_delim(forward_count_file, lazy=FALSE) -> for_data
read_delim(reverse_count_file, lazy=FALSE) -> rev_data

# Remove unwanted columns.  Calculate midpoint
for_data %>%
  select(-(`Probe Strand`:Distance)) %>%
  mutate(Position = (Start+End-1)/2) %>%
  select(-Start, -End) %>%
  select(Probe,Chromosome,Position,everything()) -> for_data

rev_data %>%
  select(-(`Probe Strand`:Distance)) %>%
  mutate(Position = (Start+End-1)/2) %>%
  select(-Start, -End) %>%
  select(Probe,Chromosome,Position,everything()) -> rev_data


# Restructure and combine
for_data %>%
  pivot_longer(cols=-(Probe:Position), names_to="Sample", values_to="For_Count") -> for_data

rev_data %>%
  pivot_longer(cols=-(Probe:Position), names_to="Sample", values_to="Rev_Count") -> rev_data

full_join(for_data,rev_data) -> data

# Clean up
rm(for_data)
rm(rev_data)

# Filtering
data %>%
  group_by(Probe) %>%
  summarise(value=mean(For_Count+Rev_Count)) %>%
  filter(value>=min_count & value <=max_count) %>%
  add_column(passed=TRUE) %>%
  select(Probe,passed) %>%
  left_join(data) %>%
  select(-passed) -> data

# Calculate directionality
data %>%
  mutate(Directionality = (Rev_Count - For_Count)/(Rev_Count + For_Count)) -> data

# Call zones
data %>%
  group_by(Sample, Chromosome) %>%
  summarise(
    call_zones(Position,Directionality,Sample[1],Chromosome[1])
  ) -> zones

zones %>%
  write_csv(paste0(output_dir,"/fork_calls.csv"))




