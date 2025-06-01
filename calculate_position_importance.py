import numpy as np
import pandas as pd
from sklearn.feature_selection import mutual_info_classif
from sklearn.metrics import classification_report, accuracy_score
from sklearn.model_selection import cross_val_score
from scipy.stats import chi2_contingency, fisher_exact
from Bio import SeqIO, AlignIO
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

class ProteinPositionAnalyzer:
    def __init__(self, gap_threshold=0.3):
        self.sequences = None
        self.labels = None
        self.gene_ids = None
        self.encoded_sequences = None
        self.results = None
        self.gap_threshold = gap_threshold
        self.original_positions = None
        self.filtered_positions = None
        self.gap_statistics = None
        
        # Amino acid to number mapping
        self.aa_to_num = {
            'A': 1, 'R': 2, 'N': 3, 'D': 4, 'C': 5,
            'Q': 6, 'E': 7, 'G': 8, 'H': 9, 'I': 10,
            'L': 11, 'K': 12, 'M': 13, 'F': 14, 'P': 15,
            'S': 16, 'T': 17, 'W': 18, 'Y': 19, 'V': 20,
            '-': 0, 'X': 0
        }
        
        # Number to amino acid mapping
        self.num_to_aa = {v: k for k, v in self.aa_to_num.items()}
        
        # Biochemical properties of amino acids
        self.aa_properties = {
            'hydrophobic': {'A', 'I', 'L', 'M', 'F', 'W', 'Y', 'V'},
            'polar': {'N', 'C', 'Q', 'S', 'T'},
            'charged': {'R', 'H', 'K', 'D', 'E'},
            'positive': {'R', 'H', 'K'},
            'negative': {'D', 'E'},
            'aromatic': {'F', 'W', 'Y'},
            'small': {'A', 'C', 'D', 'G', 'N', 'P', 'S', 'T', 'V'},
            'tiny': {'A', 'C', 'G', 'S'},
            'aliphatic': {'A', 'I', 'L', 'V'}
        }
        
    def load_labels(self, label_file):
        """Load label file"""
        labels_df = pd.read_csv(label_file, sep='\t', header=None, 
                               names=['gene_id', 'label'])
        return labels_df
    
    def load_msa(self, msa_file, format='fasta'):
        """Load multiple sequence alignment file"""
        alignment = AlignIO.read(msa_file, format)
        sequences = {}
        for record in alignment:
            sequences[record.id] = str(record.seq)
        return sequences
    
    def filter_high_gap_positions(self, sequences):
        """Filter positions with gap ratio exceeding threshold"""
        if not sequences:
            return sequences, [], []
        
        seq_length = len(sequences[0])
        n_sequences = len(sequences)
        
        gap_stats = []
        valid_positions = []
        
        print(f"Analyzing gap distribution for {seq_length} positions...")
        print(f"Gap threshold set to: {self.gap_threshold*100:.1f}%")
        
        for pos in range(seq_length):
            gap_count = sum(1 for seq in sequences if seq[pos] == '-')
            gap_ratio = gap_count / n_sequences
            
            gap_stats.append({
                'position': pos + 1,
                'gap_count': gap_count,
                'gap_ratio': gap_ratio,
                'valid': gap_ratio <= self.gap_threshold
            })
            
            if gap_ratio <= self.gap_threshold:
                valid_positions.append(pos)
        
        filtered_sequences = []
        for seq in sequences:
            filtered_seq = ''.join(seq[pos] for pos in valid_positions)
            filtered_sequences.append(filtered_seq)
        
        self.gap_statistics = pd.DataFrame(gap_stats)
        self.original_positions = list(range(1, seq_length + 1))
        self.filtered_positions = [pos + 1 for pos in valid_positions]
        
        total_positions = len(gap_stats)
        valid_positions_count = len(valid_positions)
        filtered_count = total_positions - valid_positions_count
        
        print(f"Original positions: {total_positions}")
        print(f"Valid positions: {valid_positions_count}")
        print(f"Filtered positions: {filtered_count} ({filtered_count/total_positions*100:.1f}%)")
        
        if filtered_count > 0:
            high_gap_positions = self.gap_statistics[~self.gap_statistics['valid']]['position'].tolist()
            print(f"Filtered positions: {high_gap_positions[:10]}{'...' if len(high_gap_positions) > 10 else ''}")
        
        return filtered_sequences, valid_positions, gap_stats
    
    def encode_amino_acids(self, sequences):
        """Encode amino acid sequences to numerical values"""
        encoded_seqs = []
        for seq in sequences:
            encoded_seq = [self.aa_to_num.get(aa.upper(), 0) for aa in seq]
            encoded_seqs.append(encoded_seq)
        
        return np.array(encoded_seqs)
    
    def prepare_data(self, msa_file, label_file, format='fasta'):
        """Prepare data"""
        labels_df = self.load_labels(label_file)
        sequences_dict = self.load_msa(msa_file, format)
        
        matched_sequences = []
        matched_labels = []
        matched_gene_ids = []
        
        for _, row in labels_df.iterrows():
            gene_id = row['gene_id']
            if gene_id in sequences_dict:
                matched_sequences.append(sequences_dict[gene_id])
                matched_labels.append(row['label'])
                matched_gene_ids.append(gene_id)
        
        print(f"Successfully matched {len(matched_sequences)} sequences")
        print(f"Class 0: {matched_labels.count(0)} sequences")
        print(f"Class 1: {matched_labels.count(1)} sequences")
        
        filtered_sequences, valid_positions, gap_stats = self.filter_high_gap_positions(matched_sequences)
        
        self.encoded_sequences = self.encode_amino_acids(filtered_sequences)
        self.labels = np.array(matched_labels)
        self.gene_ids = matched_gene_ids
        self.sequences = filtered_sequences
        
        return self
    
    def calculate_information_gain(self):
        """Calculate information gain"""
        mi_scores = mutual_info_classif(self.encoded_sequences, self.labels)
        return mi_scores
    
    def calculate_chi2_scores(self):
        """Calculate chi-square test scores"""
        chi2_scores = []
        p_values = []
        
        for pos in range(self.encoded_sequences.shape[1]):
            position_data = self.encoded_sequences[:, pos]
            
            unique_aas = np.unique(position_data)
            contingency_table = []
            
            for aa in unique_aas:
                class0_count = np.sum((position_data == aa) & (self.labels == 0))
                class1_count = np.sum((position_data == aa) & (self.labels == 1))
                contingency_table.append([class0_count, class1_count])
            
            contingency_table = np.array(contingency_table)
            
            if contingency_table.shape[0] > 1 and np.all(contingency_table.sum(axis=1) > 0):
                chi2_stat, p_val, _, _ = chi2_contingency(contingency_table)
                chi2_scores.append(chi2_stat)
                p_values.append(p_val)
            else:
                chi2_scores.append(0)
                p_values.append(1.0)
        
        return np.array(chi2_scores), np.array(p_values)
    
    def calculate_perfect_discrimination_score(self):
        """
        Calculate perfect discrimination score for each position
        This measures how well a position can perfectly separate the two classes
        """
        perfect_scores = []
        class_purity_scores = []
        
        for pos in range(self.encoded_sequences.shape[1]):
            position_data = self.encoded_sequences[:, pos]
            
            # Get amino acid distributions for each class
            class0_aas = position_data[self.labels == 0]
            class1_aas = position_data[self.labels == 1]
            
            # Remove gaps
            class0_aas = class0_aas[class0_aas != 0]
            class1_aas = class1_aas[class1_aas != 0]
            
            if len(class0_aas) == 0 or len(class1_aas) == 0:
                perfect_scores.append(0)
                class_purity_scores.append(0)
                continue
            
            # Calculate amino acid frequencies for each class
            class0_counts = Counter(class0_aas)
            class1_counts = Counter(class1_aas)
            
            class0_total = len(class0_aas)
            class1_total = len(class1_aas)
            
            # Method 1: Perfect discrimination score
            # Check if there are amino acids that are exclusive to one class
            class0_exclusive = set(class0_counts.keys()) - set(class1_counts.keys())
            class1_exclusive = set(class1_counts.keys()) - set(class0_counts.keys())
            
            # Calculate coverage of exclusive amino acids
            class0_exclusive_coverage = sum(class0_counts[aa] for aa in class0_exclusive) / class0_total if class0_exclusive else 0
            class1_exclusive_coverage = sum(class1_counts[aa] for aa in class1_exclusive) / class1_total if class1_exclusive else 0
            
            # Perfect discrimination score: average of exclusive coverages
            perfect_score = (class0_exclusive_coverage + class1_exclusive_coverage) / 2
            perfect_scores.append(perfect_score)
            
            # Method 2: Class purity score
            # Calculate the maximum difference in amino acid frequencies between classes
            all_aas = set(class0_counts.keys()) | set(class1_counts.keys())
            max_diff = 0
            
            for aa in all_aas:
                freq0 = class0_counts.get(aa, 0) / class0_total
                freq1 = class1_counts.get(aa, 0) / class1_total
                max_diff = max(max_diff, abs(freq0 - freq1))
            
            class_purity_scores.append(max_diff)
        
        return np.array(perfect_scores), np.array(class_purity_scores)

    def calculate_biochemical_difference_score(self):
        """Calculate biochemical property difference scores"""
        biochemical_scores = []
        
        for pos in range(self.encoded_sequences.shape[1]):
            position_data = self.encoded_sequences[:, pos]
            
            # Get amino acids for each class (excluding gaps)
            class0_aas = [self.num_to_aa[aa] for aa in position_data[self.labels == 0] if aa != 0]
            class1_aas = [self.num_to_aa[aa] for aa in position_data[self.labels == 1] if aa != 0]
            
            if len(class0_aas) == 0 or len(class1_aas) == 0:
                biochemical_scores.append(0)
                continue
            
            max_property_diff = 0
            
            # Calculate property differences
            for prop_name, prop_aas in self.aa_properties.items():
                class0_prop_freq = sum(1 for aa in class0_aas if aa in prop_aas) / len(class0_aas)
                class1_prop_freq = sum(1 for aa in class1_aas if aa in prop_aas) / len(class1_aas)
                
                prop_diff = abs(class0_prop_freq - class1_prop_freq)
                max_property_diff = max(max_property_diff, prop_diff)
            
            biochemical_scores.append(max_property_diff)
        
        return np.array(biochemical_scores)
    
    def calculate_conservation_scores(self):
        """Calculate conservation scores based on Shannon entropy"""
        conservation_scores = []
        
        for pos in range(self.encoded_sequences.shape[1]):
            position_data = self.encoded_sequences[:, pos]
            position_data = position_data[position_data != 0]
            
            if len(position_data) == 0:
                conservation_scores.append(0)
                continue
            
            aa_counts = Counter(position_data)
            total = len(position_data)
            entropy = 0
            
            for count in aa_counts.values():
                if count > 0:
                    p = count / total
                    entropy -= p * np.log2(p)
            
            max_entropy = np.log2(len(aa_counts)) if len(aa_counts) > 1 else 1
            conservation = 1 - (entropy / max_entropy) if max_entropy > 0 else 1
            conservation_scores.append(conservation)
        
        return np.array(conservation_scores)
    
    def analyze_positions(self):
        """Complete position analysis"""
        print("\nCalculating various feature importance metrics...")
        
        # Calculate various metrics
        mi_scores = self.calculate_information_gain()
        chi2_scores, p_values = self.calculate_chi2_scores()
        conservation = self.calculate_conservation_scores()
        perfect_scores, class_purity_scores = self.calculate_perfect_discrimination_score()
        biochemical_scores = self.calculate_biochemical_difference_score()
        
        # Create results DataFrame
        self.results = pd.DataFrame({
            'Original_Position': self.filtered_positions,
            'Filtered_Position': range(1, len(mi_scores) + 1),
            'Mutual_Info': mi_scores,
            'Chi2_Score': chi2_scores,
            'P_value': p_values,
            'Conservation': conservation,
            'Perfect_Discrimination': perfect_scores,
            'Class_Purity': class_purity_scores,
            'Biochemical_Difference': biochemical_scores
        })
        
        # Normalize scores
        for col in ['Mutual_Info', 'Chi2_Score', 'Perfect_Discrimination', 'Class_Purity', 'Biochemical_Difference']:
            col_data = self.results[col]
            if col_data.max() != col_data.min():
                self.results[f'{col}_norm'] = (col_data - col_data.min()) / (col_data.max() - col_data.min())
            else:
                self.results[f'{col}_norm'] = 0
        
        # Calculate combined score with higher weight for perfect discrimination
        self.results['Combined_Score'] = (
            self.results['Mutual_Info_norm'] * 0.2 + 
            self.results['Chi2_Score_norm'] * 0.2 + 
            self.results['Perfect_Discrimination_norm'] * 0.3 +  # Higher weight
            self.results['Class_Purity_norm'] * 0.15 +
            self.results['Biochemical_Difference_norm'] * 0.15
        )
        
        # Add significance markers
        self.results['Significant'] = self.results['P_value'] < 0.05
        
        # Sort by combined score
        self.results = self.results.sort_values('Combined_Score', ascending=False)
        
        return self.results
    
    def visualize_gap_statistics(self):
        """Visualize gap statistics"""
        if self.gap_statistics is None:
            print("No gap statistics available to display")
            return
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        
        # Gap ratio distribution histogram
        ax1 = axes[0, 0]
        ax1.hist(self.gap_statistics['gap_ratio'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
        ax1.axvline(x=self.gap_threshold, color='red', linestyle='--', linewidth=2, label=f'Threshold {self.gap_threshold*100:.1f}%')
        ax1.set_xlabel('Gap Ratio')
        ax1.set_ylabel('Number of Positions')
        ax1.set_title('Distribution of Gap Ratios')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Gap ratio along sequence
        ax2 = axes[0, 1]
        colors = ['red' if not valid else 'green' for valid in self.gap_statistics['valid']]
        ax2.scatter(self.gap_statistics['position'], self.gap_statistics['gap_ratio'], 
                   c=colors, alpha=0.6, s=20)
        ax2.axhline(y=self.gap_threshold, color='red', linestyle='--', linewidth=2, label=f'Threshold {self.gap_threshold*100:.1f}%')
        ax2.set_xlabel('Position')
        ax2.set_ylabel('Gap Ratio')
        ax2.set_title('Gap Ratio Distribution Along Sequence')
        ax2.legend(['Filtered', 'Valid', 'Threshold'])
        ax2.grid(True, alpha=0.3)
        
        # Valid vs filtered positions pie chart
        ax3 = axes[1, 0]
        valid_count = self.gap_statistics['valid'].sum()
        invalid_count = len(self.gap_statistics) - valid_count
        ax3.pie([valid_count, invalid_count], labels=['Valid Positions', 'Filtered Positions'], 
               autopct='%1.1f%%', colors=['lightgreen', 'lightcoral'])
        ax3.set_title('Position Filtering Statistics')
        
        # Gap count distribution
        ax4 = axes[1, 1]
        ax4.hist(self.gap_statistics['gap_count'], bins=30, alpha=0.7, color='orange', edgecolor='black')
        threshold_count = len(self.labels) * self.gap_threshold
        ax4.axvline(x=threshold_count, color='red', linestyle='--', linewidth=2, label=f'Threshold {threshold_count:.1f}')
        ax4.set_xlabel('Gap Count')
        ax4.set_ylabel('Number of Positions')
        ax4.set_title('Distribution of Gap Counts')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
        
        # Print statistics summary
        print("\nGap filtering statistics summary:")
        print("=" * 50)
        print(f"Total positions: {len(self.gap_statistics)}")
        print(f"Valid positions: {valid_count}")
        print(f"Filtered positions: {invalid_count}")
        print(f"Filtering ratio: {invalid_count/len(self.gap_statistics)*100:.1f}%")
        print(f"Average gap ratio: {self.gap_statistics['gap_ratio'].mean()*100:.2f}%")
        print(f"Maximum gap ratio: {self.gap_statistics['gap_ratio'].max()*100:.2f}%")
        
        return fig
    
    def visualize_results(self, top_n=20):
        """Visualize analysis results"""
        top_results = self.results.head(top_n)
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        
        # Score comparison
        ax1 = axes[0, 0]
        x = range(len(top_results))
        ax1.plot(x, top_results['Mutual_Info_norm'], 'o-', label='Mutual Information', alpha=0.7)
        ax1.plot(x, top_results['Chi2_Score_norm'], 's-', label='Chi-square', alpha=0.7)
        ax1.plot(x, top_results['Perfect_Discrimination_norm'], '^-', label='Perfect Discrimination', alpha=0.7)
        ax1.plot(x, top_results['Biochemical_Difference_norm'], 'v-', label='Biochemical Difference', alpha=0.7)
        ax1.plot(x, top_results['Combined_Score'], 'D-', label='Combined Score', linewidth=2)
        ax1.set_xticks(x[::2])
        ax1.set_xticklabels(top_results['Original_Position'].iloc[::2], rotation=45)
        ax1.set_xlabel('Original Position')
        ax1.set_ylabel('Normalized Score')
        ax1.legend()
        ax1.set_title('Feature Importance Scores (Top Positions)')
        ax1.grid(True, alpha=0.3)
        
        # P-value distribution
        ax2 = axes[0, 1]
        colors = ['red' if sig else 'blue' for sig in top_results['Significant']]
        ax2.bar(x, -np.log10(top_results['P_value'] + 1e-10), color=colors, alpha=0.7)
        ax2.set_xticks(x[::2])
        ax2.set_xticklabels(top_results['Original_Position'].iloc[::2], rotation=45)
        ax2.set_xlabel('Original Position')
        ax2.set_ylabel('-log10(P-value)')
        ax2.set_title('Statistical Significance')
        ax2.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='p=0.05')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Perfect discrimination vs biochemical difference
        ax3 = axes[0, 2]
        scatter = ax3.scatter(top_results['Perfect_Discrimination'], top_results['Biochemical_Difference'], 
                            c=top_results['Combined_Score'], cmap='viridis', alpha=0.7)
        ax3.set_xlabel('Perfect Discrimination Score')
        ax3.set_ylabel('Biochemical Difference Score')
        ax3.set_title('Perfect Discrimination vs Biochemical Difference')
        plt.colorbar(scatter, ax=ax3, label='Combined Score')
        ax3.grid(True, alpha=0.3)
        
        # Conservation vs discrimination
        ax4 = axes[1, 0]
        scatter = ax4.scatter(top_results['Conservation'], top_results['Combined_Score'], 
                            c=top_results['P_value'], cmap='viridis_r', alpha=0.7)
        ax4.set_xlabel('Conservation Score')
        ax4.set_ylabel('Combined Discrimination Score')
        ax4.set_title('Conservation vs Discrimination')
        plt.colorbar(scatter, ax=ax4, label='P-value')
        ax4.grid(True, alpha=0.3)
        
        # Heatmap of top positions
        ax5 = axes[1, 1]
        heatmap_data = top_results[['Mutual_Info_norm', 'Chi2_Score_norm', 'Perfect_Discrimination_norm', 
                                  'Biochemical_Difference_norm']].T
        sns.heatmap(heatmap_data, annot=False, cmap='YlOrRd', ax=ax5,
                   xticklabels=top_results['Original_Position'], 
                   yticklabels=['Mutual Info', 'Chi-square', 'Perfect Disc.', 'Biochemical Diff.'])
        ax5.set_title('Top Positions Heatmap')
        ax5.set_xlabel('Original Position')
        
        # Score distribution
        ax6 = axes[1, 2]
        score_columns = ['Perfect_Discrimination', 'Biochemical_Difference', 'Class_Purity']
        for i, col in enumerate(score_columns):
            ax6.hist(self.results[col], bins=30, alpha=0.5, label=col.replace('_', ' '))
        ax6.set_xlabel('Score Value')
        ax6.set_ylabel('Frequency')
        ax6.set_title('Distribution of Discrimination Scores')
        ax6.legend()
        ax6.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
        
        return fig
    
    def get_position_details(self, original_position):
        """Get detailed information for a specific position"""
        pos_data = self.results[self.results['Original_Position'] == original_position]
        
        if pos_data.empty:
            print(f"Position {original_position} may have been filtered or does not exist")
            if self.gap_statistics is not None:
                gap_info = self.gap_statistics[self.gap_statistics['position'] == original_position]
                if not gap_info.empty and not gap_info.iloc[0]['valid']:
                    gap_ratio = gap_info.iloc[0]['gap_ratio']
                    print(f"This position was filtered due to high gap ratio (gap ratio: {gap_ratio*100:.1f}%)")
            return None
        
        pos_result = pos_data.iloc[0]
        filtered_position = pos_result['Filtered_Position']
        pos_idx = filtered_position - 1
        
        position_data = self.encoded_sequences[:, pos_idx]
        
        class0_aas = position_data[self.labels == 0]
        class1_aas = position_data[self.labels == 1]
        
        print(f"\nDetailed analysis for position {original_position}:")
        print("=" * 50)
        print(f"Combined score: {pos_result['Combined_Score']:.4f}")
        print(f"Mutual information score: {pos_result['Mutual_Info']:.4f}")
        print(f"Chi-square score: {pos_result['Chi2_Score']:.4f}")
        print(f"Perfect discrimination score: {pos_result['Perfect_Discrimination']:.4f}")
        print(f"Biochemical difference score: {pos_result['Biochemical_Difference']:.4f}")
        print(f"P-value: {pos_result['P_value']:.2e}")
        print(f"Conservation: {pos_result['Conservation']:.4f}")
        print(f"Significant: {'Yes' if pos_result['Significant'] else 'No'}")
        
        print(f"\nAmino acid distribution in class 0:")
        class0_counts = Counter([self.num_to_aa[aa] for aa in class0_aas])
        for aa, count in sorted(class0_counts.items()):
            print(f"  {aa}: {count} ({count/len(class0_aas)*100:.1f}%)")
        
        print(f"\nAmino acid distribution in class 1:")
        class1_counts = Counter([self.num_to_aa[aa] for aa in class1_aas])
        for aa, count in sorted(class1_counts.items()):
            print(f"  {aa}: {count} ({count/len(class1_aas)*100:.1f}%)")
        
        return {
            'original_position': original_position,
            'filtered_position': filtered_position,
            'stats': pos_result,
            'class0_distribution': class0_counts,
            'class1_distribution': class1_counts
        }
    
    def save_results(self, output_file):
        """Save analysis results"""
        if self.results is not None:
            self.results.to_csv(output_file, index=False)
            print(f"Results saved to: {output_file}")
        else:
            print("No analysis results to save, please run analysis first")
    
    def save_gap_statistics(self, output_file):
        """Save gap statistics"""
        if self.gap_statistics is not None:
            self.gap_statistics.to_csv(output_file, index=False)
            print(f"Gap statistics saved to: {output_file}")
        else:
            print("No gap statistics to save")

# Usage example
def main():
    # Create analyzer with 30% gap threshold
    analyzer = ProteinPositionAnalyzer(gap_threshold=0.3)
    
    # Prepare data
    analyzer.prepare_data('test.aln.fas', 'labels.txt')
    
    # Visualize gap statistics
    analyzer.visualize_gap_statistics()
    
    # Perform analysis
    results = analyzer.analyze_positions()
    
    # Display top 30 results
    print("\nTop 30 most important positions:")
    print("=" * 120)
    print(results[['Original_Position', 'Combined_Score', 'Perfect_Discrimination', 'Biochemical_Difference',
                  'Mutual_Info', 'Chi2_Score', 'P_value', 'Significant']].head(30).to_string(index=False))
    
    # Visualize results
    analyzer.visualize_results(top_n=30)
    
    # View specific position details
    top_position = results.iloc[0]['Original_Position']
    analyzer.get_position_details(int(top_position))
    
    # Save results
    analyzer.save_results('discriminative_positions_results.csv')
    analyzer.save_gap_statistics('gap_statistics.csv')
    
    return analyzer

if __name__ == "__main__":
    analyzer = main()
