/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.filebrowsers;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

/**
 *
 * @author Juanjo Vega
 */
public class JSearchBox extends JPanel {

    protected JTextField textField = new JTextField(15);
    protected Box m_sbox = new Box(BoxLayout.X_AXIS);

    public JSearchBox(String filterText) {
        textField.setBorder(null);

        // Init UI.
        this.removeAll();
        this.setLayout(new BoxLayout(this, BoxLayout.X_AXIS));

        m_sbox.removeAll();
        m_sbox.add(Box.createHorizontalStrut(3));
        m_sbox.add(textField);
        m_sbox.add(Box.createHorizontalStrut(3));
        m_sbox.add(new CancelButton());
        m_sbox.add(Box.createHorizontalStrut(3));
        m_sbox.setBorder(BorderFactory.createLineBorder(getForeground()));
        m_sbox.setMaximumSize(new Dimension(400, 100));
        m_sbox.setPreferredSize(new Dimension(171, 20));

        Box b = new Box(BoxLayout.X_AXIS);
        b.add(new JLabel(filterText));
        b.add(Box.createHorizontalStrut(3));
        b.add(m_sbox);

        this.add(b);
    }

    public JTextField getTextField() {
        return textField;
    }

    class CancelButton extends JComponent implements MouseListener {

        private Color m_cancelColor = Color.decode("#657296"), m_cancelColorHover = Color.decode("#BEB7D8");
        private boolean hover = false;
        private int[] outline = new int[]{
            0, 0, 2, 0, 4, 2, 5, 2, 7, 0, 9, 0, 9, 2, 7, 4, 7, 5, 9, 7, 9, 9,
            7, 9, 5, 7, 4, 7, 2, 9, 0, 9, 0, 7, 2, 5, 2, 4, 0, 2, 0, 0
        };
        private int[] fill = new int[]{
            1, 1, 8, 8, 1, 2, 7, 8, 2, 1, 8, 7, 7, 1, 1, 7, 8, 2, 2, 8, 1, 8, 8, 1
        };

        public CancelButton() {
            // set button size
            Dimension d = new Dimension(10, 10);
            this.setPreferredSize(d);
            this.setMinimumSize(d);
            this.setMaximumSize(d);

            // prevent the widget from getting the keyboard focus
            this.setFocusable(false);

            // add callbacks
            this.addMouseListener(this);
        }

        @Override
        public void paintComponent(Graphics g) {
            if (hover) {    // draw fill (color depends on mouse over or not).
                g.setColor(m_cancelColorHover);
            } else {
                g.setColor(m_cancelColor);
            }
            for (int i = 0; i + 3 < fill.length; i += 4) {
                g.drawLine(fill[i], fill[i + 1], fill[i + 2], fill[i + 3]);
            }

            g.setColor(getForeground());
            for (int i = 0; i + 3 < outline.length; i += 2) {
                g.drawLine(outline[i], outline[i + 1],
                        outline[i + 2], outline[i + 3]);
            }
        }

        public void mouseClicked(MouseEvent evt) {
            textField.setText("");
        }

        public void mousePressed(MouseEvent evt) {
        }

        public void mouseReleased(MouseEvent evt) {
        }

        public void mouseEntered(MouseEvent evt) {
            hover = true;
            repaint();
        }

        public void mouseExited(MouseEvent evt) {
            hover = false;
            repaint();
        }
    }
}
