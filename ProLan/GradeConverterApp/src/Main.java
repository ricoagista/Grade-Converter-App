import converter.GradeConverter;
import converter.Result;
import java.util.Scanner;

public class Main {
    public static void main(String[] args) {
        Scanner input = new Scanner(System.in);
        GradeConverter converter = new GradeConverter();

        System.out.print("Masukkan nilai (0 - 100): ");
        double nilai = input.nextDouble();

        Result hasil = converter.convert(nilai);
        System.out.println(hasil);
    }
}
